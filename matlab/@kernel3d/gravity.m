function obj = gravity(type, gravpars, opts)
%KERNEL3D.GRAVITY   Construct a flexural-gravity surface wave kernel in 3D.
%
%   KERNEL3D.GRAVITY(type, gravpars) or KERNEL3D.GRAVITY(type, gravpars, opts),
%   where type is one of:
%
%      's'              - single layer of the free-surface Green's function,
%                         [1 1].  Which Green's function is selected by
%                         opts.green (see below).
%      'grad'           - target Cartesian gradient of the same single layer,
%                         [2 1], stacked as
%                            submat(1:2:end,:) = d/dx,
%                            submat(2:2:end,:) = d/dy.
%      'free_plate_bcs' - free-plate traces of the gravity Green's function,
%                         [2 1], stacked as
%                            row 1: M_nn + nu*M_tt   (bending moment)
%                            row 2: V_n              (Kirchhoff shear)
%                         nu is supplied in opts.nu.
%      'vol'            - volume differential operator applied to the gravity
%                         Green's function, [1 1]:
%                            (a*Laplacian^2 + b) G_S - g G_phi
%                         with [a b] supplied in opts.coefs.
%
%   GRAVPARS is the gravity parameter g (= 2*rho); the free-surface Green's
%   function has a single dispersion root rho = g/2 with residue 1, and
%   G_phi = G_S / g.
%
%   OPTS is an options struct:
%      opts.green - 'gs' (default) or 'gphi', selecting G_S or G_phi.
%                   Ignored by 'free_plate_bcs' and 'vol', which fix their
%                   own combination.
%      opts.nu    - Poisson ratio, required by 'free_plate_bcs'.
%      opts.coefs - [a b], required by 'vol'.
%
%
% See also SURFWAVE.GRAVITY.KERN, KERNEL3D.CAPILLARY, KERNEL3D.ICEFLEX

if ( nargin < 1 )
    error('KERNEL3D.GRAVITY: missing gravity kernel type.');
end
if ( nargin < 2 || isempty(gravpars) )
    error('KERNEL3D.GRAVITY: missing gravity parameter g.');
end
if ( nargin < 3 || isempty(opts) ), opts = struct(); end

green = 'gs';
if ( isfield(opts, 'green') ), green = lower(opts.green); end
if ( ~any(strcmp(green, {'gs', 'gphi'})) )
    error('KERNEL3D.GRAVITY: opts.green must be ''gs'' or ''gphi''.');
end

obj           = kernel3d();
obj.name      = 'gravity';
obj.zk        = 0;
obj.ifcomplex = 1;
obj.params.g     = gravpars;
obj.params.green = green;

% iker selectors into getnearquad_gravity_all: 0 = G_S, 1 = G_phi.
ikers = struct('gs', 0, 'gphi', 1);

switch lower(type)

    case {'s', 'single'}
        obj.type         = 's';
        obj.opdims       = [1 1];
        obj.kernel_order = -1;

        ktype = [green '_s'];
        obj.eval = @(s,t) surfwave.gravity.kern(gravpars, s, t, ktype);

        iker = ikers.(green);
        obj.getquad = @(S,eps,varargin) grav_getquad(S, eps, gravpars, ...
                          iker, 0, green, kernel3d.rsc_interleave_full(1,1), varargin{:});

    case {'grad', 'gradient'}
        obj.type         = 'grad';
        obj.opdims       = [2 1];
        obj.kernel_order = 0;

        ktype = [green 'grad_s'];
        obj.eval = @(s,t) surfwave.gravity.kern(gravpars, s, t, ktype);

        iker = ikers.(green);
        obj.getquad = @(S,eps,varargin) grav_getquad(S, eps, gravpars, ...
                          iker, 1, green, kernel3d.rsc_interleave_full(2,1), varargin{:});

    case {'free_plate_bcs', 'bcs'}
        if ( ~isfield(opts, 'nu') )
            error('KERNEL3D.GRAVITY: ''free_plate_bcs'' requires opts.nu.');
        end
        nu = opts.nu;
        obj.params.nu    = nu;
        obj.type         = 'free_plate_bcs';
        obj.opdims       = [2 1];
        obj.kernel_order = 1;
        obj.targ_fields  = {'n', 'd', 'd2'};

        obj.eval = @(s,t) surfwave.gravity.kern(gravpars, s, t, 'free_plate_bcs', nu);

        obj.getquad = [];

    case {'vol', 'volume'}
        if ( ~isfield(opts, 'coefs') )
            error('KERNEL3D.GRAVITY: ''vol'' requires opts.coefs = [a b].');
        end
        coefs = opts.coefs;
        obj.params.coefs = coefs;
        obj.type         = 'vol';
        obj.opdims       = [1 1];
        obj.kernel_order = 1;

        obj.eval = @(s,t) surfwave.gravity.kern(gravpars, s, t, 'vol', coefs);
        obj.getquad = [];

    otherwise
        error('KERNEL3D.GRAVITY: unknown gravity kernel type ''%s''.', type);

end

obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, ...
                                     obj.zk, obj.kernel_order);

end


function spmat = grav_getquad(S, eps, g, iker, ifgrad, green, ri, targinfo, opts)
%GRAV_GETQUAD  Near-quadrature correction for a gravity surface wave kernel.
%
%   Builds the RSC near pattern for (S, targinfo), calls the Fortran
%   gateway, and returns the correction as a sparse matrix.

if ( nargin < 8 || isempty(targinfo) ), targinfo = S;        end
if ( nargin < 9 || isempty(opts) ),     opts     = struct(); end

[srcvals, srccoefs, norders, ixyzs, iptype, ~] = extract_arrays(S);
npatches = S.npatches;
npts     = S.npts;

if ( isfield(opts, 'rsc') && ~isempty(opts.rsc) )
    rsc = opts.rsc;
else
    rsc = getnear(S, targinfo);
end

% On-surface targets carry patch ids; off-surface targets are flagged -1.
ntarg = size(extract_targ_array(targinfo), 2);
if ( isfield(targinfo, 'patch_id') || isprop(targinfo, 'patch_id') )
    patch_id = targinfo.patch_id;
    uvs_targ = targinfo.uvs_targ;
else
    patch_id = -ones(ntarg, 1);
    uvs_targ = zeros(2, ntarg);
end
if ( isempty(patch_id) )
    patch_id = -ones(ntarg, 1);
    uvs_targ = zeros(2, ntarg);
end

nquad = rsc.iquad(end) - 1;
% zpars(1) = g, padded to length 6
zpars = complex([g; 0; 0; 0; 0; 0]);

% scaling of the Laplace single layer added for the algebraic term
if ( strcmp(green, 'gphi') ), S3d_scal = 1; else, S3d_scal = g; end
if ( isfield(opts, 'S3d_scal') ), S3d_scal = opts.S3d_scal; end

if ( ifgrad )
    wnear = surfwave.gravity.getnearquad_gravity_grad(npatches, norders, ...
        ixyzs, iptype, npts, srccoefs, srcvals, targinfo, patch_id, ...
        uvs_targ, eps, 1, rsc.nnz, rsc.row_ptr, rsc.col_ind, rsc.iquad, ...
        rsc.rfac0, zpars, nquad, iker, S3d_scal);
else
    wnear = surfwave.gravity.getnearquad_gravity(npatches, norders, ...
        ixyzs, iptype, npts, srccoefs, srcvals, targinfo, patch_id, ...
        uvs_targ, eps, 1, rsc.nnz, rsc.row_ptr, rsc.col_ind, rsc.iquad, ...
        rsc.rfac0, zpars, nquad, iker, S3d_scal);
end

if ( size(wnear, 1) == nquad && size(wnear, 2) ~= nquad )
    wnear = wnear.';
end

spmat = conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, wnear, ri);

end

function obj = capillary(type, rts, ejs, opts)
%KERNEL3D.CAPILLARY   Construct a capillary surface wave kernel in 3D.
%
%   KERNEL3D.CAPILLARY(type, rts, ejs) or KERNEL3D.CAPILLARY(type, rts, ejs, opts).
%
%   RTS and EJS are the three dispersion roots and their partial fraction
%   residues, as returned by SURFWAVE.CAPILLARY.FIND_ROOTS_CAPILLARY.  The
%   Fortran gateway takes them packed as zpars = [rts; ejs].
%
%   Types ([1 1]):
%
%      's'       - single layer. opts.green selects G_S ('gs', default)
%                  or G_phi ('gphi').
%      'lap'     - surface Laplacian of the single layer. Near quadrature
%                  only for G_phi.
%      's3d'     - the 3D single layer applied to G_phi.
%      's3d_lap' - the Laplace 3D single layer.
%      's3d_sum' - s3d_lap + s3d.
%      'sp'      - target normal derivative of the single layer.
%      'd'       - source normal derivative of the single layer (no
%                  near quadrature).
%      'dp'      - target normal derivative of 'd' for G_S (no near
%                  quadrature).
%
%   OPTS:
%      opts.green - 'gs' (default) or 'gphi', for 's', 'lap', 'sp', 'd'.
%
% See also SURFWAVE.CAPILLARY.KERN, KERNEL3D.GRAVITY, KERNEL3D.ICEFLEX

if ( nargin < 1 )
    error('KERNEL3D.CAPILLARY: missing capillary kernel type.');
end
if ( nargin < 3 )
    error('KERNEL3D.CAPILLARY: need type, rts and ejs.');
end
if ( nargin < 4 || isempty(opts) ), opts = struct(); end

rts = rts(:);  ejs = ejs(:);
zpars = complex([rts; ejs]);

green = 'gs';
if ( isfield(opts, 'green') ), green = lower(opts.green); end
if ( ~any(strcmp(green, {'gs', 'gphi'})) )
    error('KERNEL3D.CAPILLARY: opts.green must be ''gs'' or ''gphi''.');
end

obj           = kernel3d();
obj.name      = 'capillary';
obj.ifcomplex = 1;
obj.zk        = max(abs(rts));
obj.params.rts   = rts;
obj.params.ejs   = ejs;
obj.params.green = green;

% iker selectors into getnearquad_capillary_all:
%   0 = G_S, 1 = G_phi, 3 = lap G_phi, 5 = S3d G_phi,
%   6 = Laplace S3d, 7 = 5 + 6, 8 = S'_S, 9 = S'_phi
ikers_s  = struct('gs', 0, 'gphi', 1);
ikers_sp = struct('gs', 8, 'gphi', 9);

switch lower(type)

    case {'s', 'single'}
        obj.type         = 's';
        obj.opdims       = [1 1];
        obj.kernel_order = -1;
        obj.eval    = @(s,t) surfwave.capillary.kern(rts, ejs, s, t, [green '_s']);
        obj.getquad = cap_getquad_handle(zpars, ikers_s.(green), 1, 1);

    case {'lap', 'lap_s'}
        obj.type         = 'lap';
        obj.opdims       = [1 1];
        obj.kernel_order = 1;
        obj.eval    = @(s,t) surfwave.capillary.kern(rts, ejs, s, t, ['lap_' green]);
        if ( strcmp(green, 'gphi') )
            obj.getquad = cap_getquad_handle(zpars, 3, 1, 1);
        else
            obj.getquad = [];
        end

    case {'s3d', 's3d_gphi'}
        obj.type         = 's3d';
        obj.opdims       = [1 1];
        obj.kernel_order = -1;
        obj.eval    = @(s,t) surfwave.capillary.kern(rts, ejs, s, t, 's3d_gphi');
        obj.getquad = cap_getquad_handle(zpars, 5, 1, 1);

    case {'s3d_lap'}
        obj.type         = 's3d_lap';
        obj.opdims       = [1 1];
        obj.kernel_order = -1;
        obj.eval    = @(s,t) surfwave.flex.lap3dkern(s.r, t.r);
        obj.getquad = cap_getquad_handle(zpars, 6, 1, 1);

    case {'s3d_sum', 's3d_plus_s3d_gphi'}
        obj.type         = 's3d_sum';
        obj.opdims       = [1 1];
        obj.kernel_order = -1;
        obj.eval    = @(s,t) surfwave.flex.lap3dkern(s.r, t.r) + ...
                             surfwave.capillary.kern(rts, ejs, s, t, 's3d_gphi');
        obj.getquad = cap_getquad_handle(zpars, 7, 1, 1);

    case {'d', 'double'}
        obj = cap_evalonly(obj, rts, ejs, [green '_d'], 'd', [1 1], 0);
        obj.src_fields = {'n'};

    case {'sp', 'sprime'}
        obj = cap_evalonly(obj, rts, ejs, [green '_sprime'], 'sp', [1 1], 0);
        obj.targ_fields = {'n'};
        obj.getquad = cap_getquad_handle(zpars, ikers_sp.(green), 1, 1);

    case {'dp', 'dprime'}
        obj = cap_evalonly(obj, rts, ejs, 'gs_dprime', 'dp', [1 1], 1);
        obj.src_fields = {'n'}; obj.targ_fields = {'n'};

    otherwise
        error('KERNEL3D.CAPILLARY: unknown capillary kernel type ''%s''.', type);

end

obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, ...
                                     obj.zk, obj.kernel_order);

end


function obj = cap_evalonly(obj, rts, ejs, ktype, tname, opdims, korder)
%CAP_EVALONLY  Populate an eval-only capillary type (no Fortran getquad).
obj.type         = tname;
obj.opdims       = opdims;
obj.kernel_order = korder;
obj.eval = @(s,t) surfwave.capillary.kern(rts, ejs, s, t, ktype);
obj.getquad = [];
end


function h = cap_getquad_handle(zpars, iker, m, n)
%CAP_GETQUAD_HANDLE  getquad handle for a scalar capillary Green's function.
ri = kernel3d.rsc_interleave_full(m, n);
h  = @(S,eps,varargin) cap_getquad(S, eps, zpars, iker, ri, varargin{:});
end


function spmat = cap_getquad(S, eps, zpars, iker, ri, targinfo, opts)
%CAP_GETQUAD  Near-quadrature correction for a scalar capillary kernel.
%
%   Targets are generic: targinfo may be the source surfer, another surfer,
%   or a plain struct with field .r.  Targets that lie on the source surface
%   should carry .patch_id / .uvs_targ; anything else is treated as
%   off-surface.

if ( nargin < 7 || isempty(targinfo) ), targinfo = S;        end
if ( nargin < 8 || isempty(opts) ),     opts     = struct(); end

[srcvals, srccoefs, norders, ixyzs, iptype, ~] = extract_arrays(S);

if ( isfield(opts, 'rsc') && ~isempty(opts.rsc) )
    rsc = opts.rsc;
else
    rsc = getnear(S, targinfo);
end
nquad = rsc.iquad(end) - 1;

ntarg    = size(extract_targ_array(targinfo), 2);
patch_id = [];
if ( isfield(targinfo, 'patch_id') || isprop(targinfo, 'patch_id') )
    patch_id = targinfo.patch_id;
    uvs_targ = targinfo.uvs_targ;
end
if ( isempty(patch_id) )
    patch_id = -ones(ntarg, 1);
    uvs_targ = zeros(2, ntarg);
end

wnear = surfwave.capillary.getnearquad_capillary(S.npatches, norders, ...
    ixyzs, iptype, S.npts, srccoefs, srcvals, targinfo, patch_id, ...
    uvs_targ, eps, 1, rsc.nnz, rsc.row_ptr, rsc.col_ind, rsc.iquad, ...
    rsc.rfac0, complex(zpars), nquad, iker);

if ( size(wnear, 1) == nquad && size(wnear, 2) ~= nquad )
    wnear = wnear.';
end

spmat = conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, wnear, ri);

end

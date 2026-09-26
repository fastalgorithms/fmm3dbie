function obj = flex2d(type, zpars, nu, opts)
%KERNEL3D.FLEX2D   Thin-plate (biharmonic / flexural) kernels on a surfer.
%
%   KERNEL3D.FLEX2D(type, zpars) or KERNEL3D.FLEX2D(type, zpars, nu)
%
%   The 2D thin-plate kernels of FLEX2D.KERN on a surfer in the z = 0
%   plane.
%
%   ZPARS is the plate wavenumber pair [zk1, zk2].  For the operator
%   (Delta^2 + alpha) the standard choice is zk1 = (-alpha)^(1/4),
%   zk2 = 1i*zk1.  A scalar zk is promoted to [zk, 1i*zk].
%
%   Supported types:
%
%      's'               - plate single layer, [1 1]
%
%      'free_plate_bcs'  - free-plate traces of the plate single layer,
%                          [2 1], stacked per target as
%                             row 1: M_nn + nu*M_tt  (bending moment)
%                             row 2: V_n             (Kirchhoff shear)
%                          Targets must carry n, d and d2.
%
%   NU is the Poisson ratio, required by the free-plate types.
%
%   OPTS is an options struct:
%      opts.kernel_order - override the oversampling kernel order.
%
% See also FLEX2D.KERN, FLEX2D.GETNEARQUAD, KERNEL3D.ICEFLEX, SURFERMAT

if ( nargin < 1 )
    error('KERNEL3D.FLEX2D: missing kernel type.');
end
if ( nargin < 2 || isempty(zpars) )
    error('KERNEL3D.FLEX2D: missing plate wavenumber zpars.');
end
if ( nargin < 3 ), nu = []; end
if ( nargin < 4 || isempty(opts) ), opts = struct(); end

zpars = zpars(:).';
if ( isscalar(zpars) )
    zpars = [zpars, 1i*zpars];
end
if ( numel(zpars) ~= 2 )
    error('KERNEL3D.FLEX2D: zpars must be [zk1, zk2] (or a scalar zk).');
end
zpars = complex(zpars(:));

obj              = kernel3d();
obj.name         = 'flex2d';
obj.ifcomplex    = 1;
obj.zk           = max(abs(zpars));
obj.params.zpars = zpars;
obj.params.nu    = nu;

switch lower(type)

    case {'s', 'single'}
        obj.type         = 's';
        obj.opdims       = [1 1];
        obj.kernel_order = -1;
        obj.eval         = @(s,t) flex2d.kern(zpars, s, t, 's', nu);
        ri               = kernel3d.rsc_interleave_full(1, 1);
        obj.getquad      = @(S,eps,varargin) flex2d_getquad(S, eps, ...
                               zpars, nu, 'v2v', 1, ri, varargin{:});

    case {'free_plate_bcs', 'bcs', 'v2b'}
        if ( isempty(nu) )
            error('KERNEL3D.FLEX2D: ''free_plate_bcs'' requires nu.');
        end
        obj.type         = 'free_plate_bcs';
        obj.opdims       = [2 1];
        obj.kernel_order = 1;
        obj.targ_fields  = {'n', 'd', 'd2'};
        obj.eval         = @(s,t) flex2d.kern(zpars, s, t, 'free_plate_bcs', nu);
        % 'free' returns [supp2, free2]: row 1 = moment, row 2 = shear
        ri               = kernel3d.rsc_interleave_full(2, 1);
        obj.getquad      = @(S,eps,varargin) flex2d_getquad(S, eps, ...
                               zpars, nu, 'free', 2, ri, varargin{:});

    otherwise
        error('KERNEL3D.FLEX2D: unknown plate kernel type ''%s''.', type);

end

if ( isfield(opts, 'kernel_order') ), obj.kernel_order = opts.kernel_order; end

obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, ...
                                     obj.zk, obj.kernel_order);

end


function spmat = flex2d_getquad(S, eps, zpars, nu, qtype, nker, ri, targinfo, opts)
%FLEX2D_GETQUAD  Near-quadrature correction for a plate kernel on a surfer.
%
%   Returns the accurate near-field entries (replacement values, the
%   convention SURFERMAT expects) as a sparse matrix.

if ( nargin < 8 || isempty(targinfo) ), targinfo = S;        end
if ( nargin < 9 || isempty(opts) ),     opts     = struct(); end

[srcvals, srccoefs, norders, ixyzs, iptype, ~] = extract_arrays(S);

if ( isfield(opts, 'rsc') && ~isempty(opts.rsc) )
    rsc = opts.rsc;
else
    rsc = getnear(S, targinfo);
end

nquad = rsc.iquad(end) - 1;

ntarg = size(extract_targ_array(targinfo), 2);
patch_id = [];
if ( isfield(targinfo, 'patch_id') || isprop(targinfo, 'patch_id') )
    patch_id = targinfo.patch_id;
    uvs_targ = targinfo.uvs_targ;
end
if ( isempty(patch_id) )
    patch_id = -ones(ntarg, 1);
    uvs_targ = zeros(2, ntarg);
end

dpars = nu;
if ( isempty(dpars) ), dpars = 0; end

wnear = flex2d.getnearquad(S.npatches, norders, ixyzs, iptype, S.npts, ...
    srccoefs, srcvals, targinfo, patch_id, uvs_targ, eps, 1, rsc.nnz, ...
    rsc.row_ptr, rsc.col_ind, rsc.iquad, rsc.rfac0, dpars, zpars, nquad, qtype);

wnear = reshape(wnear, nquad, nker).';

spmat = conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, wnear, ri);

end

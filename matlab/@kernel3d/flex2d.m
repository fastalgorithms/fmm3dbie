function obj = flex2d(type, zk, nu)
%KERNEL3D.FLEX2D   Construct a flexural (thin-plate) kernel in 2D.
%
%   KERNEL3D.FLEX2D(type, zk) or KERNEL3D.FLEX2D(type, zk, nu), for
%   densities on a flat surfer in the z = 0 plane, where zk is the pair
%   of plate wavenumbers [zk1, zk2] (a scalar zk is promoted to
%   [zk, 1i*zk], and zk = 0 gives the biharmonic kernel), nu is the
%   Poisson ratio, and type is:
%      's'                   - single layer G,                   [1 1]
%      'clamped_plate_bcs'   - [G; d/dn_x G],                    [2 1]
%      'supported_plate_bcs' - [G; M_nn + nu*M_tt],               [2 1]
%      'free_plate_bcs'      - [M_nn + nu*M_tt; Kirchhoff shear], [2 1]
%
%   The plate boundary condition kernels need target normals n; the
%   supported and free kernels also need target tangents d, and the free
%   kernel second derivatives d2 (or the curvature kappa). nu is required
%   for the supported and free kernels.
%
% See also FLEX2D.KERN, FLEX2D.GET_QUADRATURE_CORRECTION

if ( nargin < 2 )
    error('KERNEL3D.FLEX2D: requires type and zk arguments.');
end
if ( nargin < 3 )
    nu = [];
end

zk = zk(:).';
if ( isscalar(zk) && abs(zk) >= 1e-6 )
    zk = [zk, 1i*zk];
end

obj           = kernel3d();
obj.name      = 'flex2d';
obj.zk        = max(abs(zk));
obj.ifcomplex = 1;
obj.params.zk = zk;
obj.params.nu = nu;

switch lower(type)

    case {'s', 'single'}
        obj.type         = 's';
        obj.opdims       = [1 1];
        obj.kernel_order = -1;

    case {'clamped_plate_bcs', 'clamped'}
        obj.type         = 'clamped_plate_bcs';
        obj.opdims       = [2 1];
        obj.kernel_order = -1;
        obj.targ_fields  = {'n'};

    case {'supported_plate_bcs', 'supported'}
        obj.type         = 'supported_plate_bcs';
        obj.opdims       = [2 1];
        obj.kernel_order = -1;
        obj.targ_fields  = {'n', 'd'};

    case {'free_plate_bcs', 'free'}
        obj.type         = 'free_plate_bcs';
        obj.opdims       = [2 1];
        obj.kernel_order = -1;
        obj.targ_fields  = {'n', 'd', 'd2'};

    otherwise
        error('KERNEL3D.FLEX2D: unknown plate kernel type ''%s''.', type);

end

if ( any(strcmp(obj.type, {'supported_plate_bcs', 'free_plate_bcs'})) && isempty(nu) )
    error('KERNEL3D.FLEX2D: ''%s'' requires nu.', obj.type);
end

obj.eval    = @(s,t) flex2d.kern(zk, s, t, obj.type, nu);
nker = obj.opdims(1);
obj.getquad = @(S,eps,varargin) rsc_to_sparse( ...
    flex2d.get_quadrature_correction(S, obj.type, zk, nu, eps, varargin{:}), S, nker);
obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

% no fmm for the free plate kernel (third derivatives)
if ( strcmp(obj.type, 'free_plate_bcs') )
    obj.fmm = [];
else
    obj.fmm = @(eps,s,t,sigma) flex2d.fmm(eps, zk, s, t, obj.type, sigma, nu);
end

icheck = exist(['fmm2d.' mexext], 'file');
if ( icheck ~= 3 )
    obj.fmm = [];
end

end

function spmat = rsc_to_sparse(Q, S, nker)
%RSC_TO_SPARSE  Convert flex2d getquad RSC output to a sparse matrix.
% The plate bcs kernels store two rows per target, Q.wnear is (2,nquad).
if ( nker == 1 )
    spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear);
else
    spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear, ...
        kernel3d.rsc_interleave_full(nker, 1));
end
end

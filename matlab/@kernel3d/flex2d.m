function obj = flex2d(type, zk, nu, pfun, varargin)
%KERNEL3D.FLEX2D   Construct a flexural (thin-plate) kernel in 2D.
%
%   KERNEL3D.FLEX2D(type, zk) or KERNEL3D.FLEX2D(type, zk, nu), for
%   densities on a flat surfer in the z = 0 plane, where zk is the pair
%   of plate wavenumbers [zk1, zk2] (a scalar zk is promoted to
%   [zk, 1i*zk], and zk = 0 gives the biharmonic kernel), nu is the
%   Poisson ratio, and type is:
%      's'                   - single layer G,                   [1 1]
%      'sp'                  - sprime, d/dn_x G,                 [1 1]
%      'sg'                  - gradient, [d/dx G; d/dy G],       [2 1]
%      'clamped_plate_bcs'   - [G; d/dn_x G],                    [2 1]
%      'supported_plate_bcs' - [G; M_nn + nu*M_tt],               [2 1]
%      'free_plate_bcs'      - [M_nn + nu*M_tt; Kirchhoff shear], [2 1]
%      'varcoef'             - variable coefficient plate operator, [1 1]
%                              (alias 'flex_varalpha')
%
%   KERNEL3D.FLEX2D('varcoef', zk, nu, pfun) is alpha*I + K, where for
%   u = V[sigma] with V the 's' kernel,
%
%      Delta(alpha Delta u) - beta u - (1-nu)(alpha_xx u_yy
%         - 2 alpha_xy u_xy + alpha_yy u_xx) = alpha sigma + K[sigma].
%
%   The identity term alpha sigma is stored in obj.diag
%
%   pfun is a function handle of the target struct t (with t.r) that
%   returns a struct with fields alpha (1,nt), dalpha (2,nt), d2alpha
%   (3,nt: xx, xy, yy) and beta (1,nt). alpha and beta are not rescaled,
%   see FLEX2D.PLATE_COEFS.
%
%   The plate boundary condition kernels need target normals n; the
%   supported and free kernels also need target tangents d, and the free
%   kernel second derivatives d2 (or the curvature kappa). nu is required
%   for the supported and free kernels.
%
%   KERNEL3D.FLEX2D(..., 'ifppquad', ifppquad) returns a KERNEL3D.RADGREEN
%   kernel if true (default), otherwise uses the Fortran near quadrature.
%
% See also FLEX2D.KERN, FLEX2D.GET_QUADRATURE_CORRECTION

if ( nargin < 2 )
    error('KERNEL3D.FLEX2D: requires type and zk arguments.');
end
% positional nu and pfun, then name-value options
args = varargin;
if ( nargin >= 4 ), args = [{pfun}, args]; end
if ( nargin >= 3 ), args = [{nu}, args]; end
npos = find(cellfun(@(v) ischar(v) || isstring(v), args), 1) - 1;
if ( isempty(npos) ), npos = numel(args); end
if ( npos > 2 )
    error('KERNEL3D.FLEX2D: at most nu and pfun before the options.');
end
nu = [];  pfun = [];
if ( npos >= 1 ), nu = args{1}; end
if ( npos >= 2 ), pfun = args{2}; end
varargin = args(npos+1:end);
ifppquad = true;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'ifppquad'
            ifppquad = logical(varargin{i+1});
        otherwise
            error('KERNEL3D.FLEX2D: unknown option ''%s''.', varargin{i});
    end
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

    case {'s', 'single', 'slp'}
        obj.type         = 's';
        obj.opdims       = [1 1];
        obj.kernel_order = -1;

    case {'sp', 'sprime'}
        obj.type         = 'sp';
        obj.opdims       = [1 1];
        obj.kernel_order = -1;
        obj.targ_fields  = {'n'};

    case {'sg', 'sgrad'}
        obj.type         = 'sg';
        obj.opdims       = [2 1];
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

    case {'varcoef', 'var', 'flex_varalpha', 'varalpha'}
        obj.type         = 'varcoef';
        obj.opdims       = [1 1];
        obj.kernel_order = -1;
        if ( isempty(nu) || isempty(pfun) )
            error('KERNEL3D.FLEX2D: ''varcoef'' requires nu and pfun.');
        end
        obj.params.pfun  = pfun;
        % varcoef kernel contains its identity term
        obj.diag         = @(t) reshape(pfun(t).alpha, [], 1);

    otherwise
        error('KERNEL3D.FLEX2D: unknown plate kernel type ''%s''.', type);

end

if ( any(strcmp(obj.type, {'supported_plate_bcs', 'free_plate_bcs'})) && isempty(nu) )
    error('KERNEL3D.FLEX2D: ''%s'' requires nu.', obj.type);
end

obj.eval    = @(s,t) flex2d.kern(zk, s, t, obj.type, nu, pfun);
nker = obj.opdims(1);
obj.getquad = @(S,eps,varargin) rsc_to_sparse( ...
    getquad_rsc(S, obj.type, zk, nu, pfun, eps, varargin{:}), S, nker);
obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

% no fmm for the free plate kernel (third derivatives)
if ( strcmp(obj.type, 'free_plate_bcs') )
    obj.fmm = [];
else
    obj.fmm = @(eps,s,t,sigma) flex2d.fmm(eps, zk, s, t, obj.type, sigma, nu, pfun);
end

icheck = exist(['fmm2d.' mexext], 'file');
if ( icheck ~= 3 )
    obj.fmm = [];
end

% radgreen kernel on flex2d.green
if ( ifppquad )
    rtypes = struct('s', 's', 'sp', 'sp', 'sg', 'grad', ...
        'clamped_plate_bcs', 'clamped_bc', 'supported_plate_bcs', 'supported_bc', ...
        'free_plate_bcs', 'free_bc', 'varcoef', 'flex_varalpha');
    args = {'dim', 2, 'zk', zk};
    if ( ~isempty(nu) ), args = [args, {'nu', nu}]; end
    if ( ~isempty(pfun) ), args = [args, {'pfun', pfun}]; end
    isvar = strcmp(obj.type, 'varcoef');
    if ( ~isempty(obj.fmm) && ~isvar )
        args = [args, {'fmm', @(eps,s,pg,t,pgt) flex2d.fmmgreen(eps, zk, s, pg, t, pgt)}];
    end
    K = kernel3d(@(s,t) flex2d.green(zk, s, t), rtypes.(obj.type), args{:});
    if ( isvar )
        % third derivatives: flex2d.fmm uses Lap G_k = -k^2 G_k
        K.fmm  = obj.fmm;
        K.diag = obj.diag;
    end
    obj = K;
end

end

function Q = getquad_rsc(S, type, zk, nu, pfun, eps, targinfo, opts)
%GETQUAD_RSC  flex2d quadrature corrections, passing the plate coefficient
% function to the 'varcoef' kernel through opts.pfun.
if ( nargin < 7 ), targinfo = []; end
if ( nargin < 8 ), opts = []; end
if ( ~isempty(pfun) ), opts.pfun = pfun; end
Q = flex2d.get_quadrature_correction(S, type, zk, nu, eps, targinfo, opts);
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

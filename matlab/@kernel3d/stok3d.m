function obj = stok3d(type, coefs)
%KERNEL3D.STOK3D   Construct a Stokes kernel in 3D.
%
%   KERNEL3D.STOK3D(type) or KERNEL3D.STOK3D(type, coefs), where type is:
%      's'           - single layer (Stokeslet),    G_{ij}(x,y) = delta_ij/(8pi*r) + r_i r_j/(8pi*r^3)
%      'd','traction' - double layer (stresslet),    T_{ijk}(x,y) n_k applied to normal
%      'c'           - combined layer,              coefs(1)*S + coefs(2)*D  (default coefs=[1 1])
%
% Kernel evaluation:
%   obj.eval(srcinfo, targinfo)
%      srcinfo / targinfo are structs with fields .r (3,:) and .n (3,:).
%      Returns a (3*nt x 3*ns) matrix (opdims = [3, 3]).
%
% FMM + quadrature layer-potential eval:
%   obj.layer_eval(S, sigma, targinfo, eps)
%   obj.layer_eval(S, sigma, targinfo, eps, opts)
%      sigma is (3, ns).  Returns (3*nt, 1) column vector.
%
% Quadrature corrections:
%   obj.getquad(S, eps, varargin)
%
% Kernel orders (obj.kernel_order):  s -> -1,  d -> 0,  c -> 0
%
% See also STOK3D.KERN, STOK3D.VELOCITY.EVAL, STOK3D.VELOCITY.GET_QUADRATURE_CORRECTION

if nargin < 1
    error('KERNEL3D.STOK3D: missing Stokes kernel type.');
end

obj           = kernel3d();
obj.name      = 'stokes';
obj.opdims    = [3 3];
obj.zk        = 0;
obj.ifcomplex = 0;

switch lower(type)

    case {'s', 'single'}
        obj.type         = 's';
        obj.kernel_order = -1;

        obj.eval = @(s,t) stok_eval_reshape(stok3d.kern(s, t, 's'));

        obj.fmm = @(eps,src,targ,sigma) stok_fmm_s(eps, src, targ, sigma);

        dpars_s = [1.0; 0.0];
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            stok_layer_eval(S, sigma, targ, eps, dpars_s, varargin);
        obj.getquad  = @(S,eps,varargin) ...
                          stok_getquad(S, eps, dpars_s, varargin);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

    case {'d', 'double', 'traction'}
        obj.type         = 'd';
        obj.kernel_order = 0;

        obj.eval = @(s,t) stok_eval_reshape(stok3d.kern(s, t, 'd'));

        obj.fmm = @(eps,src,targ,sigma) stok_fmm_d(eps, src, targ, sigma);

        dpars_d = [0.0; 1.0];
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            stok_layer_eval(S, sigma, targ, eps, dpars_d, varargin);
        obj.getquad  = @(S,eps,varargin) ...
                          stok_getquad(S, eps, dpars_d, varargin);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};

    case {'c', 'combined'}
        if nargin < 2
            warning('KERNEL3D.STOK3D: missing coefs for combined layer. Defaulting to [1 1].');
            coefs = [1; 1];
        end
        coefs = coefs(:);

        obj.type         = 'c';
        obj.kernel_order = 0;
        obj.params.coefs = coefs;

        obj.eval = @(s,t) stok_eval_reshape(stok3d.kern(s, t, 'c', coefs(1), coefs(2)));

        obj.fmm = @(eps,src,targ,sigma) stok_fmm_c(eps, src, targ, sigma, coefs);

        dpars_c = coefs;
        obj.layer_eval = @(S,sigma,targ,eps,varargin) ...
                            stok_layer_eval(S, sigma, targ, eps, dpars_c, varargin);
        obj.getquad  = @(S,eps,varargin) ...
                          stok_getquad(S, eps, dpars_c, varargin);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};

    otherwise
        error('KERNEL3D.STOK3D: unknown Stokes kernel type ''%s''.', type);

end

% All Stokes velocity kernels store the symmetric 3x3 block upper triangle
% in column-major order: wnear rows = (1,1),(1,2),(1,3),(2,2),(2,3),(3,3).
obj.rsc_to_interleave = kernel3d.rsc_interleave_symmetric3();

end

% -----------------------------------------------------------------------
function mat = stok_eval_reshape(submat)
%STOK_EVAL_RESHAPE   Reshape (3,nt,3,ns) to (3*nt, 3*ns).
%
%   submat has layout (row_comp, targ, col_comp, src).
%   The target matrix M satisfies
%       M(row_comp + 3*(targ-1), col_comp + 3*(src-1)) = submat(row_comp,targ,col_comp,src)
%   which is exactly what column-major reshape gives — no permutation needed.
sz = size(submat);
if numel(sz) == 4
    nt = sz(2); ns = sz(4);
    mat = reshape(submat, [3*nt, 3*ns]);
else
    % Already scalar or flat
    mat = submat;
end
end

% -----------------------------------------------------------------------
function p = stok_layer_eval(S, sigma, targ, eps, dpars, args)
%STOK_LAYER_EVAL  Wrap stok3d.velocity.eval.
%   sigma must be (3, npts) on entry; p is returned as (3*ntarg, 1).
%   args is the varargin cell from the lambda.
if nargin < 6 || isempty(args), args = {}; end
if ~isempty(args) && isstruct(args{end})
    opts = args{end};
else
    opts = struct();
end
ntarg = size(targ.r, 2);
p_mat = stok3d.velocity.eval(S, sigma, targ, eps, dpars, opts);
p = p_mat(:);
end

% -----------------------------------------------------------------------
function Q = stok_getquad(S, eps, dpars, args)
%STOK_GETQUAD  Wrap stok3d.velocity.get_quadrature_correction.
%   args is the varargin cell from the lambda.
%   Accepts: {} (self-quadrature on S), {targinfo}, or {targinfo, opts}.
if nargin < 4 || isempty(args), args = {}; end
if length(args) >= 2
    targinfo = args{1};
    opts     = args{2};
elseif length(args) == 1
    targinfo = args{1};
    opts     = struct();
else
    targinfo = S;
    opts     = struct();
end
if isempty(targinfo)
    targinfo = S;
end
Q = stok3d.velocity.get_quadrature_correction(S, eps, dpars, targinfo, opts);
end

% -----------------------------------------------------------------------
function U = stok_fmm_s(eps, src, targ, sigma)
%STOK_FMM_S  Raw stfmm3d call for single-layer (Stokeslet) kernel.
%   sigma: (3,ns) density.  Returns (3,nt) velocity at targ.r.
srcinfo = struct();
srcinfo.sources = src.r;
srcinfo.stoklet = sigma;
U_out = stfmm3d(eps, srcinfo, 0, targ.r, 1);
U = reshape(U_out.pottarg, 3, []);
end

% -----------------------------------------------------------------------
function U = stok_fmm_d(eps, src, targ, sigma)
%STOK_FMM_D  Raw stfmm3d call for double-layer (stresslet/traction) kernel.
%   sigma: (3,ns) density; src.n: (3,ns) source normals.
%   Returns (3,nt) velocity at targ.r.
srcinfo = struct();
srcinfo.sources = src.r;
srcinfo.strslet = sigma;
srcinfo.strsvec = src.n;
U_out = stfmm3d(eps, srcinfo, 0, targ.r, 1);
U = reshape(U_out.pottarg, 3, []);
end

% -----------------------------------------------------------------------
function U = stok_fmm_c(eps, src, targ, sigma, coefs)
%STOK_FMM_C  Raw stfmm3d call for combined Stokes layer.
%   coefs(1)*S + coefs(2)*D.
%   sigma: (3,ns) density; src.n: (3,ns) source normals (used by D part).
%   Returns (3,nt) velocity at targ.r.
srcinfo = struct();
srcinfo.sources = src.r;
if coefs(1) ~= 0
    srcinfo.stoklet = coefs(1) * sigma;
end
if coefs(2) ~= 0
    srcinfo.strslet = coefs(2) * sigma;
    srcinfo.strsvec = src.n;
end
U_out = stfmm3d(eps, srcinfo, 0, targ.r, 1);
U = reshape(U_out.pottarg, 3, []);
end

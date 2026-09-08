function obj = stok3d(type, coefs)
%KERNEL3D.STOK3D   Construct a Stokes kernel in 3D.
%
%   KERNEL3D.STOK3D(type) or KERNEL3D.STOK3D(type, coefs), where type is:
%      's'            - single layer (Stokeslet),    G_{ij}(x,y) = delta_ij/(8pi*r) + r_i r_j/(8pi*r^3)
%      'd','traction'  - double layer (stresslet),    T_{ijk}(x,y) n_k(y) applied to source normal
%      'sp','sprime'   - traction of S at target,     T_{ijk}(x,y) n_k(x) (target normal)
%      'c'            - combined layer,              coefs(1)*S + coefs(2)*D  (default coefs=[1 1])
%
% Kernel evaluation:
%   obj.eval(srcinfo, targinfo)
%      srcinfo / targinfo are structs with fields .r (3,:) and .n (3,:).
%      Returns a (3*nt x 3*ns) matrix (opdims = [3, 3]).
%
% Quadrature corrections:
%   obj.getquad(S, eps, varargin)
%
% Kernel orders (obj.kernel_order):  s -> -1,  d -> 0,  sp -> 0,  c -> 0
%
% Note: 'sp'/'sprime' requires targinfo.n (target normals).
%
% See also STOK3D.KERN, STOK3D.VELOCITY.EVAL, STOK3D.VELOCITY.GET_QUADRATURE_CORRECTION,
%          STOK3D.TRACTION.EVAL, STOK3D.TRACTION.GET_QUADRATURE_CORRECTION

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

        obj.fmm = @(eps,src,targ,sigma) stok3d.fmm(eps, src, targ, 's', sigma);

        dpars_s = [1.0; 0.0];
        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          stok_getquad(S, eps, dpars_s, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

    case {'d', 'double', 'traction'}
        obj.type         = 'd';
        obj.kernel_order = 0;

        obj.eval = @(s,t) stok_eval_reshape(stok3d.kern(s, t, 'd'));

        obj.fmm = @(eps,src,targ,sigma)  stok3d.fmm(eps, src, targ, 'd', sigma);

        dpars_d = [0.0; 1.0];
        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          stok_getquad(S, eps, dpars_d, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};

    case {'sp', 'sprime', 'strac'}
        obj.type         = 'strac';
        obj.kernel_order = 0;

        obj.eval = @(s,t) stok_eval_reshape(stok3d.kern(s, t, 'sprime'));

        obj.fmm = @(eps,src,targ,sigma) stok3d.fmm(eps, src, targ, 'sp', sigma);

        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          stok_strac_getquad(S, eps, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.targ_fields = {'n'};

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

        obj.fmm = @(eps,src,targ,sigma)  stok3d.fmm(eps, src, targ, 'c', sigma, coefs);

        dpars_c = coefs;
        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          stok_getquad(S, eps, dpars_c, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};

    case {'cp', 'cprime', 'ctrac'}
        if nargin < 2
            warning('KERNEL3D.STOK3D: missing coefs for combined layer. Defaulting to [1 1].');
            coefs = [1; 1];
        end
        obj.type         = 'ctrac';
        obj.kernel_order = 1;
        dpars_c = coefs;

        obj.eval = @(s,t) stok_eval_reshape(stok3d.kern(s, t, 'cprime', coefs(1), coefs(2)));

        obj.fmm = @(eps,src,targ,sigma) stok3d.fmm(eps, src, targ, 'cp', sigma,coefs);

        obj.getquad  = @(S,eps,varargin) rsc_to_sparse( ...
                          stok_combtrac_getquad(S, eps, dpars_c, varargin{:}), S);
        obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);
        obj.src_fields = {'n'};
        obj.targ_fields = {'n'};

    otherwise
        error('KERNEL3D.STOK3D: unknown Stokes kernel type ''%s''.', type);

end

end

% -----------------------------------------------------------------------
function spmat = rsc_to_sparse(Q, S)
%RSC_TO_SPARSE  Convert Stokes getquad RSC output to a sparse matrix.
% Stokes kernels use the symmetric 3x3 upper-triangle storage.
spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear, ...
    kernel3d.rsc_interleave_symmetric3());
end

function mat = stok_eval_reshape(submat)
%STOK_EVAL_RESHAPE   Reshape (3,nt,3,ns) to (3*nt, 3*ns).
%
sz = size(submat);
if numel(sz) == 4
    nt = sz(2); ns = sz(4);
    mat = reshape(submat, [3*nt, 3*ns]);
elseif numel(sz) == 3
    nt = sz(2); ns = 1;
    mat = reshape(submat, [3*nt, 3*ns]);    
else
    % Already scalar or flat
    mat = submat;
end
end

function Q = stok_getquad(S, eps, dpars, targinfo, opts)
%STOK_GETQUAD  Wrap stok3d.velocity.get_quadrature_correction.
if nargin < 4 || isempty(targinfo), targinfo = S; end
if nargin < 5 || isempty(opts),     opts     = struct(); end
Q = stok3d.velocity.get_quadrature_correction(S, eps, dpars, targinfo, opts);
end

function Q = stok_strac_getquad(S, eps, targinfo, opts)
%STOK_SPRIME_GETQUAD  Wrap stok3d.traction.get_quadrature_correction.
if nargin < 3 || isempty(targinfo), targinfo = S; end
if nargin < 4 || isempty(opts),     opts     = struct(); end
Q = stok3d.traction.get_quadrature_correction(S, eps,[], targinfo, opts);
end

function Q = stok_combtrac_getquad(S, eps, coefs, targinfo,opts)
%STOK_SPRIME_GETQUAD  Wrap stok3d.traction.get_quadrature_correction.
if nargin < 4 || isempty(targinfo), targinfo = S; end
if nargin < 5 || isempty(opts),     opts     = struct(); end
opts.rep = 'c';
Q = stok3d.traction.get_quadrature_correction(S, eps, coefs,targinfo, opts);
end

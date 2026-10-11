function obj = eye(dvals)
%KERNEL3D.EYE   Construct an identity / "Dirac delta" kernel.
%
%       (K*sigma)(x) = D(x) * sigma(x),
%
%   K = KERNEL3D.EYE()     -> 1x1 identity delta (D = 1)
%   K = KERNEL3D.EYE(a)    -> a * delta, a numeric scalar (opdims [1 1])
%   K = KERNEL3D.EYE(A)    -> constant block delta kernel, A an m x n numeric
%                             matrix (e.g. KERNEL3D.EYE(0.5*eye(3)) for
%                             the Stokes jump)
%   K = KERNEL3D.EYE(h)    -> h(t) * delta(t-s), where h(t) returns the
%                             diagonal blocks as (m*nt x n) stacked by
%                             target, (m x n*nt) side by side, or as
%                             (m x n x nt) pages
%
%   Typical use is to add the jump term to a kernel, e.g.
%
%       kern = kernel3d('l','d') + kernel3d.eye(-0.5);
%
%   See also KERNEL3D, KERNEL3D/PLUS, KERNEL3D/TIMES, SURFERMAT.

if nargin < 1 || isempty(dvals)
    dvals = 1;
end

if isa(dvals, 'function_handle')
    D0 = probe_handle(dvals);
    m = size(D0,1); n = size(D0,2);
    dfun = @(t) coerce_diag(dvals(t), m, n, size(t.r(:,:),2));
elseif isnumeric(dvals)
    assert(ismatrix(dvals), 'KERNEL3D.EYE: matrix argument must be 2D.');
    m = size(dvals,1); n = size(dvals,2);
    A = dvals;
    dfun = @(t) repmat(A, size(t.r(:,:),2), 1);   % block A stacked per point
else
    error('KERNEL3D.EYE: argument must be a scalar, matrix, or function handle.');
end

obj           = kernel3d();
obj.name      = 'eye';
obj.type      = 'delta';
obj.opdims    = [m n];
obj.zk        = 0;
obj.ifcomplex = 0;
obj.kernel_order = -1;
obj.iszero    = false;   % zero as a smooth operator, but not as an operator
obj.src_fields  = {};
obj.targ_fields = {};
obj.diag      = dfun;

% smooth part is identically zero
obj.eval = @(s,t) zeros(m*size(t.r(:,:),2), n*size(s.r(:,:),2));

obj.fmm = @fmm_;
    function varargout = fmm_(eps, s, t, sigma) %#ok<INUSD>
        if isstruct(t), nt = size(t.r(:,:),2); else, nt = size(t(:,:),2); end
        if nargout > 0, varargout{1} = zeros(m*nt, 1); end
        if nargout > 1, varargout{2} = zeros(3, m*nt); end
        if nargout > 2
            error('KERNEL3D:eye:fmm', 'Too many output arguments.');
        end
    end

obj.getquad = @(S, eps, varargin) build_empty_quad(m, n, S, varargin);

obj.get_overs_orders = @(S,t,eps) S.norders(:);

end

function spmat = build_empty_quad(m, n, S, args)
% Return an empty (zero) sparse quadrature correction.
if ~isempty(args) && isstruct(args{1}) && isfield(args{1}, 'r')
    ntarg = size(args{1}.r(:,:), 2);
elseif ~isempty(args) && isa(args{1}, 'surfer')
    ntarg = args{1}.npts;
else
    ntarg = S.npts;
end
spmat = sparse(m*ntarg, n*S.npts);
end

function D = coerce_diag(D, m, n, nt)
% Reshape D into the (m*nt x n) stacked convention.
if m == 1 && n == 1 && numel(D) == nt
    D = D(:);                                  % scalar diag: any orientation
elseif ismatrix(D) && size(D,1) == m*nt && size(D,2) == n
    return                                     % already stacked
elseif ismatrix(D) && size(D,1) == m && size(D,2) == n*nt
    % [m x n*nt] blocks side by side -> stacked
    D = reshape(permute(reshape(D, m, n, nt), [1 3 2]), m*nt, n);
elseif size(D,1) == m && size(D,2) == n && size(D,3) == nt
    D = reshape(permute(D, [1 3 2]), m*nt, n); % [m n nt] pages -> stacked
else
    error(['KERNEL3D.EYE: diag handle must return (%d*nt x %d) stacked, ', ...
           '(%d x %d*nt) side by side, or [%d x %d x nt] pages; ', ...
           'got size [%s].'], m, n, m, n, m, n, num2str(size(D)));
end
end

function D0 = probe_handle(h)
% Evaluate h at a single random point to determine opdims.
p = probe_ptinfo();
try
    D0 = h(p);
catch
    error('KERNEL3D.EYE: unable to probe opdims from diag function handle.');
end
end

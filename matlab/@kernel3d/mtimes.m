function out = mtimes(f, g)
% * Matrix multiplication for kernel3d objects.
%
% Scalar c * K or K * c: scales eval, fmm, getquad (same as times).
%
% Left multiply M * K: M is a (p x m) matrix or function handle M(t)
%   returning (p x m x nt). Output opdims = [p, K.opdims(2)].
%
% Right multiply K * N: N is a (q x p) matrix or function handle N(s)
%   returning (q x p x ns), q = K.opdims(2). Output opdims = [K.opdims(1), p].
%
% A function handle may alternatively return the matrix forms
%   M(t): (p*nt x m) or N(s): (q x p*ns)
% or, for a pointwise scalar multiplier, its n values as a row, a column
% or a (1 x 1 x n) array.
%
% getquad returns the transformed sparse matrix (M*Q or Q*N), and the
% identity/self term diag is transformed to M(t)*D(t) or D(t)*N(t).
%
% src_fields/targ_fields are inherited from K. If M or N requires
% additional geometry fields (e.g. 'n'), add them to out.src_fields or
% out.targ_fields after calling mtimes. The handle is probed at a single
% point carrying every per-point surfer field to find its dimensions.

if isa(f, 'kernel3d') && isa(g, 'kernel3d')
    error('KERNEL3D:mtimes:invalid', ...
        'Cannot * two kernel3d objects; use + to combine.');
end

if ~isa(f, 'kernel3d')
    [f, g] = deal(g, f);
    side = 'left';
elseif ~isa(g, 'kernel3d')
    side = 'right';
else
    error('KERNEL3D:mtimes:invalid', 'Unexpected argument types.');
end

K = f;
h = g;

% scalar: delegate to times
if isnumeric(h) && isscalar(h)
    out = times(K, h);
    return;
end

% constant matrix: wrap as a constant function handle
if isnumeric(h) && ~isscalar(h)
    A = h;
    if strcmp(side, 'left')
        assert(size(A,2) == K.opdims(1), ...
            'KERNEL3D:mtimes: left matrix must have %d columns', K.opdims(1));
        p = size(A, 1);
    else
        assert(size(A,1) == K.opdims(2), ...
            'KERNEL3D:mtimes: right matrix must have %d rows', K.opdims(2));
        p = size(A, 2);
    end
    h = @(pts) repmat(A, 1, 1, size(pts.r, 2));
end

% only remaining option is a function handle
if ~isa(h, 'function_handle')
    error('KERNEL3D:mtimes:invalid', ...
        'Argument must be a scalar, matrix, or function handle.');
else
    nargfunc = nargin(h);
    assert(nargfunc==1, 'KERNEL3D:mtimes h must be a function of source or target, not both')
end

Keval    = K.eval;
Kfmm     = K.fmm;
Kgetquad = K.getquad;
Kdiag    = K.diag;
m        = K.opdims(1);
q        = K.opdims(2);

% probe h to determine output dimension p
%
% A function handle returning a (1 x 1 x n) array is treated as a
% pointwise scalar multiplier
if ~exist('p', 'var')
    try
        hval = h(probe_ptinfo());
        hr = size(hval, 1); hc = size(hval, 2);
    catch
        error('KERNEL3D:mtimes:probe', ...
            'Could not probe function handle to determine output dimension.');
    end
    if hr == 1 && hc == 1
        if strcmp(side, 'left'), p = K.opdims(1); else, p = K.opdims(2); end
    elseif strcmp(side, 'left')
        assert(hc == K.opdims(1), ...
            'KERNEL3D:mtimes: left function handle must return matrices with %d columns', K.opdims(1));
        p = hr;
    else
        assert(hr == K.opdims(2), ...
            'KERNEL3D:mtimes: right function handle must return matrices with %d rows', K.opdims(2));
        p = hc;
    end
    h = normalize_handle(h, hr, hc, side);
end

if K.iszero
    if strcmp(side, 'left'), out = kernel3d.zeros([p, q]);
    else,                    out = kernel3d.zeros([m, p]); end
    return;
end

out      = K;
out.type = ['custom_', K.type];
out.name = ['custom ', K.name];

if strcmp(side, 'left')
    out.opdims  = [p, q];
    out.eval    = @eval_left;
    out.fmm     = set_if_exist(Kfmm,     @fmm_left);
    out.getquad = set_if_exist(Kgetquad, @getquad_left);
    out.diag    = set_if_exist(Kdiag,    @diag_left);
else
    out.opdims  = [m, p];
    out.eval    = @eval_right;
    out.fmm     = set_if_exist(Kfmm,     @fmm_right);
    out.getquad = set_if_exist(Kgetquad, @getquad_right);
    out.diag    = set_if_exist(Kdiag,    @diag_right);
end

    function out = apply_left(fval, X)
        if size(fval,1) == 1 && size(fval,2) == 1
            out = fval .* X;
        else
            out = pagemtimes(fval, X);
        end
    end

    function out = apply_right(X, fval)
        if size(fval,1) == 1 && size(fval,2) == 1
            out = X .* fval;
        else
            out = pagemtimes(X, fval);
        end
    end

% left-multiply: h(t) * K(s,t)

    function vals = eval_left(s, t)
        nt   = size(t.r, 2);
        ns   = size(s.r, 2);
        fval = h(t);
        Kmat = Keval(s, t);
        K3   = permute(reshape(Kmat, m, nt, q*ns), [1 3 2]);
        out3 = apply_left(fval, K3);
        vals = reshape(permute(out3, [1 3 2]), p*nt, q*ns);
    end

    function out = fmm_left(eps, s, t, sigma)
        nt    = size(t.r, 2);
        fval  = h(t);
        inner = Kfmm(eps, s, t, sigma);
        out   = reshape(apply_left(fval, reshape(inner, m, 1, nt)), p*nt, 1);
    end

    function Q = getquad_left(S, eps, varargin)
        Qinner = Kgetquad(S, eps, varargin{:});
        if ~isempty(varargin) && isstruct(varargin{1})
            targ = varargin{1};
        else
            targ = S;
        end
        nt   = size(targ.r, 2);
        fval = h(targ);
        Q3   = permute(reshape(full(Qinner), m, nt, q*S.npts), [1 3 2]);
        Q    = sparse(reshape(permute(apply_left(fval, Q3), [1 3 2]), p*nt, q*S.npts));
    end

    function D = diag_left(t)
        % stacked (m*nt x q) diag -> pages, left-multiply, restack
        nt = size(t.r, 2);
        D3 = permute(reshape(Kdiag(t), m, nt, q), [1 3 2]);
        D  = reshape(permute(apply_left(h(t), D3), [1 3 2]), p*nt, q);
    end

% right-multiply: K(s,t) * h(s)

    function vals = eval_right(s, t)
        ns   = size(s.r, 2);
        nt   = size(t.r, 2);
        fval = h(s);
        Kmat = Keval(s, t);
        K3   = reshape(Kmat, m*nt, q, ns);
        vals = reshape(apply_right(K3, fval), m*nt, p*ns);
    end

    function out = fmm_right(eps, s, t, sigma)
        ns     = size(s.r, 2);
        fval   = h(s);
        sig_in = reshape(apply_left(fval, reshape(sigma, p, 1, ns)), q, ns);
        out    = Kfmm(eps, s, t, sig_in);
    end

    function Q = getquad_right(S, eps, varargin)
        Qinner = Kgetquad(S, eps, varargin{:});
        ns   = S.npts;
        fval = h(S);
        if ~isempty(varargin) && isstruct(varargin{1})
            nt = size(varargin{1}.r, 2);
        else
            nt = ns;
        end
        Q3 = reshape(full(Qinner), m*nt, q, ns);
        Q  = sparse(reshape(apply_right(Q3, fval), m*nt, p*ns));
    end

    function D = diag_right(t)
        % the delta term has s = t, so h is evaluated at the targets
        nt = size(t.r, 2);
        D3 = permute(reshape(Kdiag(t), m, nt, q), [1 3 2]);
        D  = reshape(permute(apply_right(D3, h(t)), [1 3 2]), m*nt, p);
    end

end

function hn = normalize_handle(h, hr, hc, side)
% Function handles may return either a tensor (hr x hc x n) or
% a stacked matrix:
%   left,  h(t): (hr*n x hc), stacked by target
%   right, h(s): (hr x hc*n), stacked by source
% Wrap h so that it always returns the tensor form.
    function fval = hwrap(pts)
        fval = h(pts);
        n = size(pts.r, 2);
        if hr == 1 && hc == 1 && numel(fval) == n
            % pointwise scalar: row, column or pages
            fval = reshape(fval, 1, 1, n);
            return
        end
        if n <= 1 || size(fval, 3) ~= 1
            return
        end
        if strcmp(side, 'left') && size(fval,1) == hr*n && size(fval,2) == hc
            fval = permute(reshape(fval, hr, n, hc), [1 3 2]);
        elseif strcmp(side, 'right') && size(fval,1) == hr && size(fval,2) == hc*n
            fval = reshape(fval, hr, hc, n);
        elseif size(fval,1) == hr && size(fval,2) == hc
            % constant (hr x hc) returned: broadcast over points
            fval = repmat(fval, 1, 1, n);
        end
    end
hn = @hwrap;
end

function out = set_if_exist(cond, val)
if isa(cond, 'function_handle')
    out = val;
else
    out = [];
end
end

function f = times(f, g)
% .* Pointwise (elementwise) multiplication for kernel3d class
%
% K.*A and A.*K (equivalent) where A is numeric and K is a kernel3d with
% opdims [m q]. The product is taken elementwise on each source-target
% pair. A may be
%   - a scalar,
%   - a column vector (m x 1),
%   - a row vector (1 x q),
%   - an m x q matrix.
%
% Singleton dimensions of either A or K.opdims are expanded, just like the
% usual .* for matrices.
%
% eval, getquad and diag are always adjusted. The fmm is adjusted when A
% is a scalar, a column vector or a row vector; it is set to [] if A is a
% general matrix.

if ~isa(f, 'kernel3d')
    f = times(g, f);
    return
end

if isa(g, 'kernel3d')
    error('KERNEL3D:times:invalid', 'Cannot .* two kernel3d objects');
end

if ~isnumeric(g)
    error('KERNEL3D:times:invalid', ...
        'F or G must be numeric and the other a kernel3d class object');
end

if ~ismatrix(g)
    error('KERNEL3D:times:invalid', 'numeric factor must be a 2D array');
end

if isscalar(g)
    f = times_scalar(f, g);
    return
end

A = g;
m = f.opdims(1); q = f.opdims(2);
[a1, a2] = size(A);
if ~((a1 == m || a1 == 1 || m == 1) && (a2 == q || a2 == 1 || q == 1))
    error('KERNEL3D:times:dims', ...
        'KERNEL3D:times: size of A (%d x %d) incompatible with kernel opdims (%d x %d)', ...
        a1, a2, m, q);
end
P = max(a1, m); Q = max(a2, q);

if f.iszero || all(A(:) == 0)
    f = kernel3d.zeros([P, Q]);
    return
end

Keval    = f.eval;
Kfmm     = f.fmm;
Kgetquad = f.getquad;
Kdiag    = f.diag;
A4 = reshape(A, a1, 1, a2, 1);

    function vals = apply_mat(Kmat)
        % Kmat is (m*nt) x (q*n); apply A blockwise with expansion
        nt = size(Kmat,1)/m;
        n  = size(Kmat,2)/q;
        K4 = reshape(full(Kmat), m, nt, q, n);
        vals = reshape(A4 .* K4, P*nt, Q*n);
    end

    function vals = apply_sparse(Kmat)
        % sparse version of apply_mat, only touches the nonzeros
        nt = size(Kmat,1)/m;
        n  = size(Kmat,2)/q;
        [i, j, v] = find(Kmat);
        it = floor((i-1)/m); ik = i - it*m;
        jt = floor((j-1)/q); jk = j - jt*q;
        I = []; J = []; V = [];
        for pp = 1:P
            if m == 1, sel_r = true(size(ik)); else, sel_r = (ik == pp); end
            for qq = 1:Q
                if q == 1, sel = sel_r; else, sel = sel_r & (jk == qq); end
                aval = A(min(pp, a1), min(qq, a2));
                I = [I; it(sel)*P + pp]; %#ok<AGROW>
                J = [J; jt(sel)*Q + qq]; %#ok<AGROW>
                V = [V; aval * v(sel)];  %#ok<AGROW>
            end
        end
        vals = sparse(I, J, V, P*nt, Q*n);
    end

    function vals = eval_(s, t)
        vals = apply_mat(Keval(s, t));
    end

    function out = fmm_col(eps, s, t, sigma)
        % A is a column vector: post-multiply the fmm output
        u  = Kfmm(eps, s, t, sigma);
        nt = numel(u)/m;
        out = reshape(A(:) .* reshape(u, m, nt), P*nt, 1);
    end

    function out = fmm_row(eps, s, t, sigma)
        % A is a row vector: pre-multiply the density
        ns  = numel(sigma)/Q;
        sig = A(:) .* reshape(sigma, Q, ns);
        if q == 1 && Q > 1
            % K has a single input channel: [a1 K, a2 K, ...] sigma
            sig = sum(sig, 1);
        end
        out = Kfmm(eps, s, t, reshape(sig, q*ns, 1));
    end

    function Qm = getquad_(S, eps, varargin)
        Qm = apply_sparse(sparse(Kgetquad(S, eps, varargin{:})));
    end

if isa(Keval, 'function_handle')
    f.eval = @eval_;
else
    f.eval = [];
end
if isa(Kfmm, 'function_handle') && a2 == 1
    f.fmm = @fmm_col;
elseif isa(Kfmm, 'function_handle') && a1 == 1
    f.fmm = @fmm_row;
else
    f.fmm = [];
end
if isa(Kgetquad, 'function_handle')
    f.getquad = @getquad_;
else
    f.getquad = [];
end

if isa(Kdiag, 'function_handle')
    f.diag = @diag_;
end

    function D = diag_(t)
        % apply A to each (m x q) point block of the stacked diag
        Dm = Kdiag(t);
        nt = size(Dm,1)/m;
        D3 = A .* permute(reshape(Dm, m, nt, q), [1 3 2]);
        D  = reshape(permute(D3, [1 3 2]), P*nt, Q);
    end

f.opdims = [P, Q];
f.type = ['custom_', f.type];
f.name = ['custom ', f.name];

end

function f = times_scalar(f, g)
if f.iszero || g == 0
    f = kernel3d.zeros(f.opdims);
    return
end

if isa(f.eval, 'function_handle')
    Keval = f.eval;
    f.eval = @(varargin) g * Keval(varargin{:});
else
    f.eval = [];
end

if isa(f.fmm, 'function_handle')
    Kfmm = f.fmm;
    f.fmm = @(varargin) g * Kfmm(varargin{:});
else
    f.fmm = [];
end

if isa(f.getquad, 'function_handle')
    Kgetquad = f.getquad;
    f.getquad = @(S, eps, varargin) kernel3d.scalequad( ...
        Kgetquad(S, eps, varargin{:}), g);
else
    f.getquad = [];
end

if isa(f.diag, 'function_handle')
    Kdiag = f.diag;
    f.diag = @(t) g * Kdiag(t);
end

if isnan(g)
    f.isnan = true;
end
end

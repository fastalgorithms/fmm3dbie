function out = times(f, g)
% .* Multiplication for kernel3d objects.
%
% Scalar:   c .* K  or  K .* c
%   Scales eval, fmm, getquad by c.
%
% Left (target) multiply:   f .* K
%   f(t) returns (p x m x nt), or a 2D array of shape
%   (nt x p*m) or (p*m x nt)
%   Output opdims = [p, K.opdims(2)].
%
% Right (source) multiply:  K .* f
%   f(s) returns (K.opdims(2) x p x ns), or a 2D array of shape
%   (ns x K.opdims(2)*p) or (K.opdims(2)*p x ns).
%   Output opdims = [K.opdims(1), p].
%
% Constant matrix:  A .* K  or  K .* A  where A is numeric (not scalar)
%   Wrapped as a constant function handle.
%
% getquad returns the outer sparse matrix (F*Q_inner or Q_inner*F).
%
% Note: src_fields/targ_fields are inherited from K. If f requires
% additional geometry fields (e.g. 'n'), add them to out.src_fields or
% out.targ_fields.

if isa(f, 'kernel3d') && isa(g, 'kernel3d')
    error('KERNEL3D:times:invalid', ...
        'Cannot .* two kernel3d objects; use + to combine.');
end

if ~isa(f, 'kernel3d')
    [f, g] = deal(g, f);
    side = 'left';
elseif ~isa(g, 'kernel3d')
    side = 'right';
else
    error('KERNEL3D:times:invalid', 'Unexpected argument types.');
end

K = f;
h = g;

% Scalar
if isnumeric(h) && isscalar(h)
    out = K;
    if isa(out.eval, 'function_handle')
        out.eval = @(varargin) h * out.eval(varargin{:});
    end
    if isa(out.fmm, 'function_handle')
        out.fmm = @(varargin) h * out.fmm(varargin{:});
    else
        out.fmm = [];
    end
    if isa(out.getquad, 'function_handle')
        Kgetquad = out.getquad;
        out.getquad = @(S, eps, varargin) kernel3d.scalequad( ...
            Kgetquad(S, eps, varargin{:}), h);
    else
        out.getquad = [];
    end
    return;
end

% Constant matrix: wrap as a constant function handle.
if isnumeric(h) && ~isscalar(h)
    A = h;
    if strcmp(side, 'left')
        assert(size(A,2) == K.opdims(1), ...
            'KERNEL3D:times: left matrix must have %d columns', K.opdims(1));
        p = size(A, 1);
    else
        assert(size(A,1) == K.opdims(2), ...
            'KERNEL3D:times: right matrix must have %d rows', K.opdims(2));
        p = size(A, 2);
    end
    h = @(pts) repmat(A, 1, 1, size(pts.r, 2));
end

% Function handle
if ~isa(h, 'function_handle')
    error('KERNEL3D:times:invalid', ...
        'Argument must be a scalar, matrix, or function handle.');
else
    nargfunc = nargin(h);
    assert(nargfunc==1, 'KERNEL3D:times h must be a function of source or target, not both')
end

Keval    = K.eval;
Kfmm     = K.fmm;
Kgetquad = K.getquad;
m        = K.opdims(1);
q        = K.opdims(2);

% Probe h to determine output dimension p (skipped if already set above).
if ~exist('p', 'var')
    try
        probe.r  = randn(3,1); probe.n  = randn(3,1);
        probe.du = randn(3,1); probe.dv = randn(3,1);
        hval = h(probe);
        rc = size(hval,1)*size(hval,2)*size(hval,3);
        if strcmp(side, 'left')
            mq = m;
        else
            mq = q;
        end
        if ndims(hval) == 3 || rc == 1
            % genuine 3D tensor, or an unambiguous scalar
            if strcmp(side, 'left')
                p = size(hval, 1);
            else
                p = size(hval, 2);
            end
        else
            % 2D output
            p = rc / mq;
        end
    catch
        error('KERNEL3D:times:probe', ...
            'Could not probe function handle to determine output dimension.');
    end
end

out      = K;
out.type = ['custom_', K.type];
out.name = ['custom ', K.name];

if strcmp(side, 'left')
    out.opdims  = [p, q];
    out.eval    = @eval_left;
    out.fmm     = set_if_exist(Kfmm,     @fmm_left);
    out.getquad = set_if_exist(Kgetquad, @getquad_left);
else
    out.opdims  = [m, p];
    out.eval    = @eval_right;
    out.fmm     = set_if_exist(Kfmm,     @fmm_right);
    out.getquad = set_if_exist(Kgetquad, @getquad_right);
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

    function fval3 = to3d(fval, rows, cols, npts)
        % Handle h's that output a (rows x cols x npts) tensor,
        % or as npts x (rows*cols) or (rows*cols) x npts matrices
        if ndims(fval) == 3 && size(fval,3) == npts
            fval3 = fval;
            return
        end
        assert(ismatrix(fval), ...
            'KERNEL3D:times: h output must be 2D or 3D.');
        rc = rows*cols;
        if isequal(size(fval), [npts, rc])
            fval3 = permute(reshape(fval, npts, rows, cols), [2 3 1]);
        elseif isequal(size(fval), [rc, npts])
            fval3 = reshape(fval, rows, cols, npts);
        elseif npts == 1 && numel(fval) == rc
            fval3 = reshape(fval, rows, cols, 1);
        else
            error('KERNEL3D:times:shape', ...
                ['h output has shape %s, expected (%d x %d x %d), ', ...
                 '(%d x %d), or (%d x %d)'], ...
                mat2str(size(fval)), rows, cols, npts, npts, rc, rc, npts);
        end
    end

% Left-multiply  (h(t) post-multiplies the output)

    function vals = eval_left(s, t)
        nt   = size(t.r, 2);
        ns   = size(s.r, 2);
        fval = to3d(h(t), p, m, nt);                % (p x m x nt)
        Kmat = Keval(s, t);                        % (m*nt x q*ns)
        K3   = permute(reshape(Kmat, m, nt, q*ns), [1 3 2]);  % (m x q*ns x nt)
        out3 = apply_left(fval, K3);               % (p x q*ns x nt)
        vals = reshape(permute(out3, [1 3 2]), p*nt, q*ns);
    end

    function out = fmm_left(eps, s, t, sigma)
        nt    = size(t.r, 2);
        fval  = to3d(h(t), p, m, nt);               % (p x m x nt)
        inner = Kfmm(eps, s, t, sigma);            % (m*nt x 1)
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
        fval = to3d(h(targ), p, m, nt);             % (p x m x nt)
        Q3   = permute(reshape(full(Qinner), m, nt, q*S.npts), [1 3 2]);
        Q    = sparse(reshape(permute(apply_left(fval, Q3), [1 3 2]), p*nt, q*S.npts));
    end

% Right-multiply  (h(s) pre-multiplies the density)

    function vals = eval_right(s, t)
        ns   = size(s.r, 2);
        nt   = size(t.r, 2);
        fval = to3d(h(s), q, p, ns);                % (q x p x ns)
        Kmat = Keval(s, t);                        % (m*nt x q*ns)
        K3   = reshape(Kmat, m*nt, q, ns);
        vals = reshape(apply_right(K3, fval), m*nt, p*ns);
    end

    function out = fmm_right(eps, s, t, sigma)
        ns     = size(s.r, 2);
        fval   = to3d(h(s), q, p, ns);              % (q x p x ns)
        sig_in = reshape(apply_left(fval, reshape(sigma, p, 1, ns)), q, ns);
        out    = Kfmm(eps, s, t, sig_in);
    end

    function Q = getquad_right(S, eps, varargin)
        Qinner = Kgetquad(S, eps, varargin{:});
        ns   = S.npts;
        fval = to3d(h(S), q, p, ns);                % (q x p x ns)
        if ~isempty(varargin) && isstruct(varargin{1})
            nt = size(varargin{1}.r, 2);
        else
            nt = ns;
        end
        Q3 = reshape(full(Qinner), m*nt, q, ns);
        Q  = sparse(reshape(apply_right(Q3, fval), m*nt, p*ns));
    end

end

% --------------------------------------------------------------------------
% Helper
% --------------------------------------------------------------------------

function out = set_if_exist(cond, val)
if isa(cond, 'function_handle')
    out = val;
else
    out = [];
end
end

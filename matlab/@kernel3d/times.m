function out = times(f, g)
% .* Multiplication for kernel3d objects.
%
% Scalar:   c .* K  or  K .* c
%   Scales eval, fmm, layer_eval, getquad by c.
%
% Left (target) multiply:   f .* K
%   f(t) returns (p x m x nt) or (1 x 1 x nt) for pointwise scaling.
%   Output opdims = [p, K.opdims(2)].
%
% Right (source) multiply:  K .* f
%   f(s) returns (q x p x ns) or (1 x 1 x ns), q = K.opdims(2).
%   Output opdims = [K.opdims(1), p].
%
% Constant matrix:  A .* K  or  K .* A  where A is numeric (not scalar)
%   Wrapped as a constant function handle.
%
% getquad returns the outer sparse matrix (F*Q_inner or Q_inner*F).
%
% layer_eval precomp protocol:
%   - Sparse precomp flows DOWN through the times-wrapper chain.
%   - params.is_times_wrapper = true marks a kernel as a times-wrapper so
%     inner wrappers know to keep passing sparse rather than converting to RSC.
%
% Note: src_fields/targ_fields are inherited from K. If f requires
% additional geometry fields (e.g. 'n'), add them to out.src_fields or
% out.targ_fields after calling times.

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
    if isa(out.layer_eval, 'function_handle')
        out.layer_eval = @(varargin) h * out.layer_eval(varargin{:});
    else
        out.layer_eval = [];
    end
    if isa(out.getquad, 'function_handle')
        Kgetquad = out.getquad;
        Kri      = out.rsc_to_interleave;
        out.getquad = @(S, eps, varargin) kernel3d.scalequad( ...
            Kgetquad(S, eps, varargin{:}), S, h, Kri);
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
    assert(nargfunc==2, 'KERNEL3D:times h must be a function of source or target, not both')

end

Keval       = K.eval;
Kfmm        = K.fmm;
Klayer_eval = K.layer_eval;
Kgetquad    = K.getquad;
Kri         = K.rsc_to_interleave;   % inner kernel's rsc descriptor
m           = K.opdims(1);
q           = K.opdims(2);

% Probe h to determine output dimension p (skipped if already set above).
if ~exist('p', 'var')
    try
        probe.r  = randn(3,1); probe.n  = randn(3,1);
        probe.du = randn(3,1); probe.dv = randn(3,1);
        hval = h(probe);
        if strcmp(side, 'left')
            p = size(hval, 1);
        else
            p = size(hval, 2);
        end
    catch
        error('KERNEL3D:times:probe', ...
            'Could not probe function handle to determine output dimension.');
    end
end

out         = K;
out.type    = ['custom_', K.type];
out.name    = ['custom ', K.name];
% Mark this kernel as a times-wrapper so that nested layer_eval calls
% know to pass sparse precomp further down.
out.params.is_times_wrapper = true;

if strcmp(side, 'left')
    out.opdims            = [p, q];
    out.rsc_to_interleave = kernel3d.rsc_interleave_full(p, q);
    out.eval              = @eval_left;
    out.fmm               = set_if_exist(Kfmm,       @fmm_left);
    out.layer_eval        = set_if_exist(Klayer_eval, @layer_eval_left);
    out.getquad           = set_if_exist(Kgetquad,    @getquad_left);
else
    out.opdims            = [m, p];
    out.rsc_to_interleave = kernel3d.rsc_interleave_full(m, p);
    out.eval              = @eval_right;
    out.fmm               = set_if_exist(Kfmm,       @fmm_right);
    out.layer_eval        = set_if_exist(Klayer_eval, @layer_eval_right);
    out.getquad           = set_if_exist(Kgetquad,    @getquad_right);
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

    function Fpinv = compute_pinv(fval)
        if size(fval,1) == 1 && size(fval,2) == 1
            Fpinv = 1 ./ fval;
        else
            Fpinv = pagepinv(fval);
        end
    end

% Left-multiply  (h(t) post-multiplies the output)

    function vals = eval_left(s, t)
        nt   = size(t.r, 2);
        ns   = size(s.r, 2);
        fval = h(t);                               % (p x m x nt)
        Kmat = Keval(s, t);                        % (m*nt x q*ns)
        K3   = permute(reshape(Kmat, m, nt, q*ns), [1 3 2]);  % (m x q*ns x nt)
        out3 = apply_left(fval, K3);               % (p x q*ns x nt)
        vals = reshape(permute(out3, [1 3 2]), p*nt, q*ns);
    end

    function out = fmm_left(eps, s, t, sigma)
        nt    = size(t.r, 2);
        fval  = h(t);                              % (p x m x nt)
        inner = Kfmm(eps, s, t, sigma);            % (m*nt x 1)
        out   = reshape(apply_left(fval, reshape(inner, m, 1, nt)), p*nt, 1);
    end

    function out = layer_eval_left(S, sigma, targ, eps, varargin)
        nt   = size(targ.r, 2);
        fval = h(targ);                            % (p x m x nt)

        if has_sparse_precomp(varargin)
            % Q_outer = F_left * Q_inner (sparse, p*nt x q*ns).
            % Strip F_left by pagepinv to recover Q_inner (sparse, m*nt x q*ns).
            opts    = varargin{end};
            Q_outer = opts.precomp_quadrature;
            ns      = S.npts;
            Fpinv   = compute_pinv(fval);          % (m x p x nt)
            Q3      = permute(reshape(full(Q_outer), p, nt, q*ns), [1 3 2]);
            Qin3    = apply_left(Fpinv, Q3);       % (m x q*ns x nt)
            Q_inner = sparse(reshape(permute(Qin3, [1 3 2]), m*nt, q*ns));

            if inner_is_times_wrapper()
                % Inner layer_eval is also a times-wrapper: pass sparse down.
                opts.precomp_quadrature = Q_inner;
            else
                % Terminal: convert to RSC for the base Fortran layer_eval.
                rsc_near = getnear(S, targ);
                Q_rsc    = conv_spmat_to_rsc(S, Q_inner, Kri);
                Q_rsc.targinfo     = targ;
                Q_rsc.rfac         = rsc_near.rfac;
                Q_rsc.wavenumber   = K.zk;
                Q_rsc.kernel_order = K.kernel_order;
                Q_rsc.format       = 'rsc';
                opts.precomp_quadrature = Q_rsc;
            end
            varargin{end} = opts;
        end

        inner = Klayer_eval(S, sigma, targ, eps, varargin{:});
        out   = reshape(apply_left(fval, reshape(inner(:), m, 1, nt)), p*nt, 1);
    end

    function Q = getquad_left(S, eps, varargin)
        Qinner = to_sparse(Kgetquad(S, eps, varargin{:}), S, Kri);
        if ~isempty(varargin) && isstruct(varargin{1})
            targ = varargin{1};
        else
            targ = S;
        end
        nt   = size(targ.r, 2);
        fval = h(targ);                            % (p x m x nt)
        Q3   = permute(reshape(full(Qinner), m, nt, q*S.npts), [1 3 2]);
        Q    = sparse(reshape(permute(apply_left(fval, Q3), [1 3 2]), p*nt, q*S.npts));
    end

% Right-multiply  (h(s) pre-multiplies the density)

    function vals = eval_right(s, t)
        ns   = size(s.r, 2);
        nt   = size(t.r, 2);
        fval = h(s);                               % (q x p x ns)
        Kmat = Keval(s, t);                        % (m*nt x q*ns)
        K3   = reshape(Kmat, m*nt, q, ns);
        vals = reshape(apply_right(K3, fval), m*nt, p*ns);
    end

    function out = fmm_right(eps, s, t, sigma)
        ns     = size(s.r, 2);
        fval   = h(s);                             % (q x p x ns)
        sig_in = reshape(apply_left(fval, reshape(sigma, p, 1, ns)), q, ns);
        out    = Kfmm(eps, s, t, sig_in);
    end

    function out = layer_eval_right(S, sigma, targ, eps, varargin)
        ns   = size(S.r, 2);
        fval = h(S);                               % (q x p x ns)
        sig_in = reshape(apply_left(fval, reshape(sigma, p, 1, ns)), q, ns);

        if has_sparse_precomp(varargin)
            % Q_outer = Q_inner * F_right (sparse, m*nt x p*ns).
            % Strip F_right by pagepinv to recover Q_inner (sparse, m*nt x q*ns).
            opts    = varargin{end};
            Q_outer = opts.precomp_quadrature;
            nt      = size(targ.r, 2);
            Fpinv   = compute_pinv(fval);          % (p x q x ns)
            Q3      = reshape(full(Q_outer), m*nt, p, ns);
            Qin3    = apply_right(Q3, Fpinv);      % (m*nt x q x ns)
            Q_inner = sparse(reshape(Qin3, m*nt, q*ns));

            if inner_is_times_wrapper()
                % Inner layer_eval is also a times-wrapper: pass sparse down.
                opts.precomp_quadrature = Q_inner;
            else
                % Terminal: convert to RSC for the base Fortran layer_eval.
                rsc_near = getnear(S, targ);
                Q_rsc    = conv_spmat_to_rsc(S, Q_inner, Kri);
                Q_rsc.targinfo     = targ;
                Q_rsc.rfac         = rsc_near.rfac;
                Q_rsc.wavenumber   = K.zk;
                Q_rsc.kernel_order = K.kernel_order;
                Q_rsc.format       = 'rsc';
                opts.precomp_quadrature = Q_rsc;
            end
            varargin{end} = opts;
        end

        out = Klayer_eval(S, sig_in, targ, eps, varargin{:});
    end

    function tf = inner_is_times_wrapper()
        % True if Klayer_eval belongs to a times-wrapped kernel (not a base kernel).
        tf = isfield(K.params, 'is_times_wrapper') && K.params.is_times_wrapper;
    end

    function Q = getquad_right(S, eps, varargin)
        Qinner = to_sparse(Kgetquad(S, eps, varargin{:}), S, Kri);
        ns   = S.npts;
        fval = h(S);                               % (q x p x ns)
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
% Helpers
% --------------------------------------------------------------------------

function Qsp = to_sparse(Q, S, ri)
if isstruct(Q) && isfield(Q, 'row_ptr')
    Qsp = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear, ri);
else
    Qsp = Q;
end
end

function tf = has_sparse_precomp(varargin_cell)
tf = ~isempty(varargin_cell) && isstruct(varargin_cell{end}) && ...
     isfield(varargin_cell{end}, 'precomp_quadrature') && ...
     issparse(varargin_cell{end}.precomp_quadrature);
end

function out = set_if_exist(cond, val)
if isa(cond, 'function_handle')
    out = val;
else
    out = [];
end
end

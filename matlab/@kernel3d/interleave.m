function K = interleave(kerns)
%KERNEL3D.INTERLEAVE   Create a kernel3d from a matrix of kernels by interleaving.
%
%   K = kernel3d.interleave(kerns)
%
%   kerns is an (m x n) matrix of kernel3d objects. The resulting kernel K
%   has opdims = [sum(rowdims), sum(coldims)] and combines the sub-kernels
%   so that:
%     - K.eval(s,t) returns the (opdims(1)*nt) x (opdims(2)*ns) matrix with
%       each (k,l) block placed at the interleaved row/column indices.
%     - K.fmm, K.layer_eval, K.getquad are wired analogously.
%
%   See also KERNEL3D.PLUS, KERNEL3D.TIMES

[m, n] = size(kerns);

rowdims   = arrayfun(@(k) kerns(k,1).opdims(1), 1:m);
coldims   = arrayfun(@(l) kerns(1,l).opdims(2), 1:n);
rowstarts = [0 cumsum(rowdims)];
colstarts = [0 cumsum(coldims)];
opdims    = [sum(rowdims) sum(coldims)];

% Worst-case kernel order across all sub-kernels
kernel_order = max(arrayfun(@(k) kerns(k).kernel_order, 1:numel(kerns)));

% -------------------------------------------------------------------------
%  Helper: build interleaved row/col index cells for given nt, ns
% -------------------------------------------------------------------------
    function [ridx, cidx] = make_idx(nt, ns)
        ridx = cell(m, 1);
        cidx = cell(n, 1);
        for k = 1:m
            ridx{k} = (rowstarts(k)+1:rowstarts(k+1)).' + (0:nt-1)*opdims(1);
            ridx{k} = ridx{k}(:).';
        end
        for l = 1:n
            cidx{l} = (colstarts(l)+1:colstarts(l+1)).' + (0:ns-1)*opdims(2);
            cidx{l} = cidx{l}(:).';
        end
    end

% -------------------------------------------------------------------------
%  eval
% -------------------------------------------------------------------------
    function out = eval_(s, t)
        [~, ns] = size(s.r);
        [~, nt] = size(t.r);
        [ridx, cidx] = make_idx(nt, ns);
        out = zeros(opdims(1)*nt, opdims(2)*ns);
        for k = 1:m
            for l = 1:n
                out(ridx{k}, cidx{l}) = kerns(k,l).eval(s, t);
            end
        end
    end

% -------------------------------------------------------------------------
%  fmm
% -------------------------------------------------------------------------
    function varargout = fmm_(eps, s, t, sigma)
        if isstruct(t)
            [~, nt] = size(t.r);
        else
            [~, nt] = size(t);
        end
        [~, ns] = size(s.r);
        [ridx, cidx] = make_idx(nt, ns);

        % Call each sub-kernel FMM and accumulate
        pot  = zeros(opdims(1)*nt, 1);
        if nargout > 1
            grad = zeros(3, opdims(1)*nt);
        end

        for k = 1:m
            for l = 1:n
                if nargout > 1
                    [p, g] = kerns(k,l).fmm(eps, s, t, sigma(cidx{l}));
                    grad(:, ridx{k}) = grad(:, ridx{k}) + g;
                else
                    p = kerns(k,l).fmm(eps, s, t, sigma(cidx{l}));
                end
                pot(ridx{k}) = pot(ridx{k}) + p(:);
            end
        end

        varargout{1} = pot;
        if nargout > 1
            varargout{2} = grad;
        end
    end

% -------------------------------------------------------------------------
%  layer_eval  -- signature: layer_eval(S, sigma, targinfo, eps, ...)
% -------------------------------------------------------------------------
    function out = layer_eval_(S, sigma, targinfo, eps, varargin)
        if isstruct(targinfo)
            nt = size(targinfo.r, 2);
        else
            nt = size(targinfo, 2);
        end
        ns = size(sigma, 2);

        [ridx, cidx] = make_idx(nt, ns);

        out = zeros(opdims(1)*nt, 1);
        for k = 1:m
            for l = 1:n
                sig_l = sigma(cidx{l});
                p = kerns(k,l).layer_eval(S, sig_l, targinfo, eps, varargin{:});
                out(ridx{k}) = out(ridx{k}) + p(:);
            end
        end
    end

% -------------------------------------------------------------------------
%  getquad  -- signature: getquad(S, eps, targinfo, ...)
%  Returns a sparse matrix (sum of all sub-kernel quad corrections).
% -------------------------------------------------------------------------
    function Q = getquad_(S, eps, varargin)
        Q = sparse(0);
        for k = 1:m
            for l = 1:n
                Qkl_raw = kerns(k,l).getquad(S, eps, varargin{:});
                % Convert rsc struct to sparse if needed
                if isstruct(Qkl_raw) && isfield(Qkl_raw, 'row_ptr')
                    Qkl = conv_rsc_to_spmat(S, Qkl_raw.row_ptr, Qkl_raw.col_ind, Qkl_raw.wnear);
                else
                    Qkl = Qkl_raw;
                end
                % Qkl is (opdims_kl(1)*nt x opdims_kl(2)*ns); embed it into
                % the full (opdims(1)*nt x opdims(2)*ns) sparse matrix.
                % We do this by building the permutation from the sub-block
                % indices.  Since nt and ns come from the quad correction
                % size, infer them.
                [nrows_kl, ncols_kl] = size(Qkl);
                nt_kl = nrows_kl / rowdims(k);
                ns_kl = ncols_kl / coldims(l);
                [ridx_kl, cidx_kl] = make_idx(nt_kl, ns_kl);
                % Embed into full sparse matrix
                [ii, jj, vv] = find(Qkl);
                % Map local row/col indices back to interleaved global indices
                ii_global = ridx_kl{k}(ii).';
                jj_global = cidx_kl{l}(jj).';
                if isempty(Q) || all(size(Q) == [1 1])
                    Q = sparse(ii_global, jj_global, vv, ...
                        opdims(1)*nt_kl, opdims(2)*ns_kl);
                else
                    Q = Q + sparse(ii_global, jj_global, vv, ...
                        size(Q,1), size(Q,2));
                end
            end
        end
    end

% -------------------------------------------------------------------------
%  get_overs_orders  -- elementwise max over all sub-kernels
% -------------------------------------------------------------------------
    function novers = get_overs_orders_(S, t, eps)
        novers = [];
        for k = 1:numel(kerns)
            nov_k = kerns(k).get_overs_orders(S, t, eps);
            if isempty(novers)
                novers = nov_k;
            else
                novers = max(novers, nov_k);
            end
        end
    end

% -------------------------------------------------------------------------
%  Union of src_fields / targ_fields across all sub-kernels
% -------------------------------------------------------------------------
all_src_fields  = {};
all_targ_fields = {};
for ki = 1:numel(kerns)
    all_src_fields  = union(all_src_fields,  kerns(ki).src_fields);
    all_targ_fields = union(all_targ_fields, kerns(ki).targ_fields);
end
if isempty(all_src_fields),  all_src_fields  = []; end
if isempty(all_targ_fields), all_targ_fields = []; end

% -------------------------------------------------------------------------
%  Assemble output kernel
% -------------------------------------------------------------------------
K           = kernel3d();
K.name      = 'custom';
K.type      = 'interleaved';
K.opdims    = opdims;
K.kernel_order = kernel_order;
K.ifcomplex   = max(arrayfun(@(k) kerns(k).ifcomplex, 1:numel(kerns)));
K.src_fields  = all_src_fields;
K.targ_fields = all_targ_fields;

K.eval = @eval_;

if all(arrayfun(@(k) isa(kerns(k).fmm, 'function_handle'), 1:numel(kerns)))
    K.fmm = @fmm_;
else
    K.fmm = [];
end

if all(arrayfun(@(k) isa(kerns(k).layer_eval, 'function_handle'), 1:numel(kerns)))
    K.layer_eval = @layer_eval_;
else
    K.layer_eval = [];
end

if all(arrayfun(@(k) isa(kerns(k).getquad, 'function_handle'), 1:numel(kerns)))
    K.getquad = @getquad_;
else
    K.getquad = [];
end

if all(arrayfun(@(k) isa(kerns(k).get_overs_orders, 'function_handle'), 1:numel(kerns)))
    K.get_overs_orders = @get_overs_orders_;
else
    K.get_overs_orders = [];
end

end

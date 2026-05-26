function obj = tangent_kern(kernin, ids, jds)
%KERNEL3D.TANGENT_KERN   Wrap a kernel so vector blocks are expressed in
%   the local tangent basis (du, dv).
%
%   OBJ = KERNEL3D.TANGENT_KERN(KERNIN, IDS, JDS)
%
%   IDS - 3 x ki matrix; each column gives the indices of a source vector
%         block. Unlisted indices are scalar and pass through unchanged.
%   JDS - 3 x ko matrix; same for target.
%
%   Let Ps     = tangent_block(s, ids, ..., true)   % (src_tan*ns) x (src_cart*ns), contraction
%       Ps_exp = tangent_block(s, ids, ..., false)  % (src_cart*ns) x (src_tan*ns), expansion
%       Pt     = tangent_block(t, jds, ..., true)   % (tgt_tan*nt) x (tgt_cart*nt), contraction
%       Pt_exp = tangent_block(t, jds, ..., false)  % (tgt_cart*nt) x (tgt_tan*nt), expansion
%
%   eval:       Pt * K_cart(s,t) * Ps_exp
%   fmm:        Pt * K_fmm(s,t, Ps_exp * sigma)
%   layer_eval: Pt * K_le(S, Ps_exp * sigma, targ, eps, ...)
%   getquad:    Pt * Q_cart * Ps_exp
%
%   For layer_eval, if opts.precomp_quadrature is present it was produced
%   by the wrapped getquad (tangent basis). It is mapped back to Cartesian
%   via Pt_exp * Q_tan * Ps, converted to rsc, and passed to the inner
%   layer_eval.
%
%   OBJ.opdims = [tgt_opdim - ko, src_opdim - ki].
%
%   See also TANGENT_BLOCK

ki = size(ids, 2);
ko = size(jds, 2);

src_opdim_in  = kernin.opdims(2);
tgt_opdim_in  = kernin.opdims(1);
src_opdim_out = src_opdim_in - ki;
tgt_opdim_out = tgt_opdim_in - ko;

src_scalar_rows = setdiff((1:src_opdim_in)', ids(:));
tgt_scalar_rows = setdiff((1:tgt_opdim_in)', jds(:));

    function out = eval_(s, t)
        Ps_exp = tangent_block(s, ids, src_scalar_rows, src_opdim_out, false);
        Pt     = tangent_block(t, jds, tgt_scalar_rows, tgt_opdim_out, true);
        out = Pt * kernin.eval(s, t) * Ps_exp;
    end

    function out = fmm_(eps, s, t, sigma)
        ns     = size(s.r, 2);
        nt     = size(t.r, 2);
        Ps_exp = tangent_block(s, ids, src_scalar_rows, src_opdim_out, false);
        Pt     = tangent_block(t, jds, tgt_scalar_rows, tgt_opdim_out, true);
        sigma_cart = reshape(Ps_exp * sigma, src_opdim_in, ns);
        out = Pt * reshape(kernin.fmm(eps, s, t, sigma_cart), tgt_opdim_in*nt, 1);
    end

    function out = layer_eval_(S, sigma, targ, eps, varargin)
        Ps_exp = tangent_block(S,    ids, src_scalar_rows, src_opdim_out, false);
        Pt     = tangent_block(targ, jds, tgt_scalar_rows, tgt_opdim_out, true);

        if ~isempty(varargin) && isstruct(varargin{end}) && ...
                isfield(varargin{end}, 'precomp_quadrature')
            opts  = varargin{end};
            Q_tan = opts.precomp_quadrature;
            if isstruct(Q_tan), Q_tan = conv_rsc_to_spmat(S,Q_tan.row_ptr,Q_tan.col_ind,Q_tan.wnear, kernin.rsc_to_interleave); end
            Ps     = tangent_block(S,    ids, src_scalar_rows, src_opdim_out, true);
            Pt_exp = tangent_block(targ, jds, tgt_scalar_rows, tgt_opdim_out, false);
            Q_cart = conv_spmat_to_rsc(S,Pt_exp * Q_tan * Ps,kernin.rsc_to_interleave);
            
            % TODO we could get rfac faster than this
            rsc_near = getnear(S, targ);

            Q_cart.targinfo     = targ;
            Q_cart.format       = 'rsc';  Q_cart.wavenumber = kernin.zk;
            Q_cart.kernel_order = kernin.kernel_order; Q_cart.rfac    = rsc_near.rfac;
            opts.precomp_quadrature = Q_cart;
            varargin{end} = opts;
        end

        ns_S = size(S.r, 2);
        nt   = size(targ.r, 2);
        sigma_cart = reshape(Ps_exp * sigma(:), src_opdim_in, ns_S);
        out = Pt * reshape(kernin.layer_eval(S, sigma_cart, targ, eps, varargin{:}), tgt_opdim_in*nt, 1);
    end

    function Q = getquad_(S, eps, varargin)
        Q_cart = kernin.getquad(S, eps, varargin{:});
        if ~issparse(Q_cart)
            Q_cart = conv_rsc_to_spmat(S, Q_cart.row_ptr, Q_cart.col_ind, Q_cart.wnear,kernin.rsc_to_interleave);
        end
        if ~isempty(varargin) && isstruct(varargin{1})
            targ = varargin{1};
        else
            targ = S;
        end
        Ps_exp = tangent_block(S,    ids, src_scalar_rows, src_opdim_out, false);
        Pt     = tangent_block(targ, jds, tgt_scalar_rows, tgt_opdim_out, true);
        Q = Pt * Q_cart * Ps_exp;
    end

obj = kernel3d();
obj.name         = kernin.name;
obj.type         = [kernin.type, '_tangent'];
obj.opdims       = [tgt_opdim_out, src_opdim_out];
obj.zk           = kernin.zk;
obj.ifcomplex    = kernin.ifcomplex;
obj.kernel_order = kernin.kernel_order;
obj.src_fields   = union(kernin.src_fields,  {'du', 'dv'});
obj.targ_fields  = union(kernin.targ_fields, {'du', 'dv'});
obj.eval         = @eval_;

if isa(kernin.fmm, 'function_handle')
    obj.fmm = @fmm_;
end
if isa(kernin.layer_eval, 'function_handle')
    obj.layer_eval = @layer_eval_;
end
if isa(kernin.getquad, 'function_handle')
    obj.getquad = @getquad_;
end
if isa(kernin.get_overs_orders, 'function_handle')
    obj.get_overs_orders = kernin.get_overs_orders;
end

end

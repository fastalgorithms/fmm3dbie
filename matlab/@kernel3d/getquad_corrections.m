function A = getquad_corrections(K, S, eps, targinfo, nover)
%KERNEL3D.GETQUAD_CORRECTIONS   Near quadrature minus the oversampled
%   smooth rule, for kernel K on surface S at targets targinfo.
%
%   A = kernel3d.getquad_corrections(K, S, eps, targinfo, nover)
%
%   nover is the oversampling order of the smooth rule (NaN for none).

if ~ismethod(S, 'oversample') || any(isnan(nover))
    nover = S.norders;
end

if K.getquad_sub
    A = K.getquad(S, eps, targinfo, struct('corrections', true, 'nover', nover));
    return
end

A = K.getquad(S, eps, targinfo);
rsc = getnear(S, targinfo);
wsmth = smooth_sparse_quad(K, targinfo, S, rsc.row_ptr, rsc.col_ind, nover);
A = A - conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, wsmth, ...
    kernel3d.rsc_interleave_full(K.opdims(1), K.opdims(2)));

end

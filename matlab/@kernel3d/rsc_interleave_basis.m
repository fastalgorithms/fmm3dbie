function ri = rsc_interleave_basis(m, n, nker, entries)
%RSC_INTERLEAVE_BASIS  rsc_to_interleave descriptor for basis-expansion kernels.
%
%   ri = kernel3d.rsc_interleave_basis(m, n, nker, entries)
%
%   For kernels where wnear(nker, nquad) stores scalar basis functions that
%   combine with coefficients to fill a (m x n) block:
%
%     B(row_id, col_id) += coef * wnear(ker_id, :)
%
%   entries is a struct array, each element with fields:
%     .ker_id  - which wnear row (1-based)
%     .row_id  - block row (1-based, in 1..m)
%     .col_id  - block col (1-based, in 1..n)
%     .coef    - scalar coefficient (may be complex)

ri.type    = 'basis';
ri.m       = m;
ri.n       = n;
ri.nker    = nker;
ri.entries = entries;

end

function ri = rsc_interleave_full(m, n)
%RSC_INTERLEAVE_FULL  rsc_to_interleave descriptor for full m-by-n block kernels.
%
%   ri = kernel3d.rsc_interleave_full(m, n)
%
%   For kernels with opdims = [m, n]: wnear has nker = m*n rows stored in
%   column-major order within the block, i.e.
%     wnear(k, :)  <->  block entry (row, col)
%   where k = row + m*(col-1),  row in 1..m,  col in 1..n.
%
%   Returns a struct with fields:
%     .type    = 'full'
%     .nker    = m*n
%     .row_ids = (m*n, 1) row index within block for each wnear row (1-based)
%     .col_ids = (m*n, 1) col index within block for each wnear row (1-based)

nker = m * n;
% Column-major ordering: row index varies fastest.
% wnear(k) <-> block entry (mod(k-1, m)+1,  ceil(k/m))
[rows, cols] = meshgrid(1:m, 1:n);   % rows(j,i)=i, cols(j,i)=j; (:) gives col-major
ri.type    = 'full';
ri.nker    = nker;
ri.row_ids = rows(:);   % (m*n, 1)  — row index for wnear row k
ri.col_ids = cols(:);   % (m*n, 1)  — col index for wnear row k

end

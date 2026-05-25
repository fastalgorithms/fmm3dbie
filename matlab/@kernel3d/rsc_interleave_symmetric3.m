function ri = rsc_interleave_symmetric3()
%RSC_INTERLEAVE_SYMMETRIC3  rsc_to_interleave descriptor for symmetric 3x3 block kernels.
%
%   ri = kernel3d.rsc_interleave_symmetric3()
%
%   For kernels with opdims = [3, 3] whose quadrature stores only the upper
%   triangle of the symmetric block in column-major order (6 entries):
%
%     wnear(1,:) <-> (1,1)
%     wnear(2,:) <-> (1,2)
%     wnear(3,:) <-> (1,3)
%     wnear(4,:) <-> (2,2)
%     wnear(5,:) <-> (2,3)
%     wnear(6,:) <-> (3,3)
%
%   This matches the Stokes velocity kernel layout (stok_comb_vel.f90).
%
%   Returns a struct with fields:
%     .type         = 'symmetric'
%     .nker         = 6
%     .row_ids      = (6,1) row index within block for each wnear row
%     .col_ids      = (6,1) col index within block for each wnear row
%     .sym_row_ids  = (3,1) row indices of the off-diagonal reflected entries
%     .sym_col_ids  = (3,1) col indices of the off-diagonal reflected entries
%     .sym_ker_ids  = (3,1) which wnear row each reflected entry mirrors

% Upper triangle, column-major order (col varies slowest):
% (1,1),(1,2),(1,3),(2,2),(2,3),(3,3)
ri.type    = 'symmetric';
ri.nker    = 6;
ri.row_ids = [1; 1; 1; 2; 2; 3];
ri.col_ids = [1; 2; 3; 2; 3; 3];

% Off-diagonal reflected entries (lower triangle): (2,1),(3,1),(3,2)
% Each mirrors wnear row 2, 3, 5 respectively
ri.sym_row_ids = [2; 3; 3];
ri.sym_col_ids = [1; 1; 2];
ri.sym_ker_ids = [2; 3; 5];

end

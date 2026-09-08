function ri = rsc_interleave_symmetric3()
%RSC_INTERLEAVE_SYMMETRIC3  rsc_to_interleave descriptor for symmetric 3x3 block kernels.
%
%   ri = kernel3d.rsc_interleave_symmetric3()
%
%   For kernels with opdims = [3, 3] whose quadrature stores only the upper
%   triangle of the symmetric block in column-major order (6 entries):
%
%     wnear(1,:) <-> (1,1)
%     wnear(2,:) <-> (1,2) and reflected (2,1)
%     wnear(3,:) <-> (1,3) and reflected (3,1)
%     wnear(4,:) <-> (2,2)
%     wnear(5,:) <-> (2,3) and reflected (3,2)
%     wnear(6,:) <-> (3,3)
%
%   This matches the Stokes velocity kernel layout (stok_comb_vel.f90).

% Upper triangle entries (column-major): (1,1),(1,2),(1,3),(2,2),(2,3),(3,3)
upper_rows = [1; 1; 1; 2; 2; 3];
upper_cols = [1; 2; 3; 2; 3; 3];

% Lower triangle reflections of the three off-diagonal entries:
%   wnear(2) -> (2,1),  wnear(3) -> (3,1),  wnear(5) -> (3,2)
refl_ker  = [2; 3; 5];
refl_rows = [2; 3; 3];
refl_cols = [1; 1; 2];

nker = 6;
ne   = nker + numel(refl_ker);   % 9 total entries
e(ne) = struct('ker_id',[], 'row_id',[], 'col_id',[], 'coef',[]);
for k = 1:nker
    e(k) = struct('ker_id', k, 'row_id', upper_rows(k), 'col_id', upper_cols(k), 'coef', 1);
end
for s = 1:numel(refl_ker)
    e(nker+s) = struct('ker_id', refl_ker(s), 'row_id', refl_rows(s), 'col_id', refl_cols(s), 'coef', 1);
end

ri = kernel3d.rsc_interleave_full(3, 3, nker, e);

end

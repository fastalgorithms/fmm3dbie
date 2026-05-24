function ri = rsc_interleave_scalar()
%RSC_INTERLEAVE_SCALAR  rsc_to_interleave descriptor for scalar kernels.
%
%   ri = kernel3d.rsc_interleave_scalar()
%
%   For kernels with opdims = [1,1]: wnear has nker=1 row, and each entry
%   maps directly to the single (1,1) block position.
%
%   Returns a struct with fields:
%     .type     = 'scalar'
%     .nker     = 1
%     .row_inds = 1
%     .col_inds = 1

ri.type     = 'scalar';
ri.nker     = 1;
ri.row_inds = 1;
ri.col_inds = 1;

end

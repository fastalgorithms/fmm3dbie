function Q = addquad(Qf, Qg, S, sign, ri)
%KERNEL3D.ADDQUAD   Add (or subtract) two quadrature correction structs.
%
%   Q = kernel3d.addquad(Qf, Qg, S)          computes Qf + Qg
%   Q = kernel3d.addquad(Qf, Qg, S, -1)      computes Qf - Qg
%   Q = kernel3d.addquad(Qf, Qg, S, sign, ri) uses ri as the
%       rsc_to_interleave descriptor when converting rsc structs to sparse.
%
%   Each input Q may be either an rsc struct (with fields row_ptr, col_ind,
%   wnear) or a sparse matrix.  Both are converted to sparse matrices,
%   combined, and returned as a sparse matrix.

if nargin < 4 || isempty(sign)
    sign = 1;
end
if nargin < 5
    ri = [];
end

% Convert Qf to sparse
if isstruct(Qf) && isfield(Qf, 'row_ptr')
    Sf = conv_rsc_to_spmat(S, Qf.row_ptr, Qf.col_ind, Qf.wnear, ri);
else
    Sf = Qf;
end

% Convert Qg to sparse
if isstruct(Qg) && isfield(Qg, 'row_ptr')
    Sg = conv_rsc_to_spmat(S, Qg.row_ptr, Qg.col_ind, Qg.wnear, ri);
else
    Sg = Qg;
end

Q = Sf + sign * Sg;
end

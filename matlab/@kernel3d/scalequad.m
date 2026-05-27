function Q = scalequad(Qf, S, c, ri)
%KERNEL3D.SCALEQUAD   Scale a quadrature correction struct by a scalar.
%
%   Q = kernel3d.scalequad(Qf, S, c)      returns c * Qf as a sparse matrix.
%   Q = kernel3d.scalequad(Qf, S, c, ri)  uses ri as the rsc_to_interleave
%       descriptor when converting an rsc struct to sparse.
%
%   Qf may be either an rsc struct (with fields row_ptr, col_ind, wnear)
%   or a sparse matrix.

if nargin < 4
    ri = [];
end

% Convert to sparse
if isstruct(Qf) && isfield(Qf, 'row_ptr')
    Sf = conv_rsc_to_spmat(S, Qf.row_ptr, Qf.col_ind, Qf.wnear, ri);
else
    Sf = Qf;
end

Q = c * Sf;
end

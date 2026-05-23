function Q = scalequad(Qf, S, c)
%KERNEL3D.SCALEQUAD   Scale a quadrature correction struct by a scalar.
%
%   Q = kernel3d.scalequad(Qf, S, c)   returns c * Qf as a sparse matrix.
%
%   Qf may be either an rsc struct (with fields row_ptr, col_ind, wnear)
%   or a sparse matrix.

% Convert to sparse
if isstruct(Qf) && isfield(Qf, 'row_ptr')
    Sf = conv_rsc_to_spmat(S, Qf.row_ptr, Qf.col_ind, Qf.wnear);
else
    Sf = Qf;
end

Q = c * Sf;
end

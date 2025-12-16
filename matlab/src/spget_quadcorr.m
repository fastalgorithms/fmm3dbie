
function A = spget_quadcorr(I, J, P, S)
% SPGET_quadcorr(I,J,P,S) computes entries of a sparse matrix of near-field
% corrections that should be added to the kernel matrix, as used in
% AFUN. P is the permutation required for extracting near
% quadrature correction, and S stores the sparse matrix corresponding
% to the quadrature correction
m_ = length(I);
n_ = length(J);
[I_sort,E] = sort(I);
P(I_sort) = E;
A = zeros(m_,n_);
[I,J,S_] = find(S(:,J));
I = cast(I, "int64");
I_sort = cast(I_sort, "int64");
idx = ismembc(I,I_sort);
I = I(idx);
J = J(idx);
S_ = S_(idx);
A(P(I) + (J - 1)*m_) = S_;
end
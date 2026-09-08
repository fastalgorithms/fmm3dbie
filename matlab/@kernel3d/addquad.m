function Q = addquad(Qa, Qb, c)
%KERNEL3D.ADDQUAD   Add or subtract two sparse quadrature correction matrices.
%
%   Q = kernel3d.addquad(Qa, Qb)      computes Qa + Qb
%   Q = kernel3d.addquad(Qa, Qb, -1)  computes Qa - Qb

if nargin < 3 || isempty(c)
    c = 1;
end
Q = Qa + c * Qb;
end

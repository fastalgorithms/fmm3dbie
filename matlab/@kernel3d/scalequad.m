function Q = scalequad(Qf, c)
%KERNEL3D.SCALEQUAD   Scale a sparse quadrature correction matrix by a scalar.
%
%   Q = kernel3d.scalequad(Qf, c)  returns c * Qf.

Q = c * Qf;
end

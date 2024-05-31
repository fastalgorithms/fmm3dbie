%  This file contains detailed info about the individual integral
%  representations for solving the perfect electric conductor
%  boundary value problem for Maxwell equations
%
%  The PDE takes the form
%  1v) \nabla \times E =  ik H
%  2v) \nabla \cdot  E =     0
%  3v) \nabla \times H = -ik E
%  4v) \nabla \cdot  H =     0
%  
%  where E is the electric field, H is the magnetic field, 
%  and k is the wavenumber
%
%  The PEC boundary conditions are given by
%  1b) n \times (E + E_in) = 0
%  2b) n \cdot  (E + E_in) = \rho
%  3b) n \times (H + H_in) = J
%  4b) n \cdot  (H + H_in) = 0
%
%  where (E_in, H_in) are the incoming electric and magnetic 
%  fields, \rho is the surface current, and J is the 
%  surface current.
%
%
%  NRCCIE (Non-resonant charge-current integral equation)
%  
%  Representation:
%  H = \nabla \times S_{k}[J]
%  E = ik S_{k}[J] - \nabla S_{k} [\rho]
%
%  Integral equations solve for [J, \rho] and obtaine dby imposing
%
%  n \times (H + H_in) - \alpha (n \times n \times (E + E_in)) = J
%  n \cdot (E + E_in) + \alpha/ik \nabla \cdot E               = \rho
%
%
%  and are given by
%
%  J/2 - n \times \nabla S_{k} J + 
%  \alpha (n \times n \times (ik S_{k} J - \nabla S_{k}[\rho])) = 
%   n \times H_in - \alpha n \times n \times E_in
%
%  \rho/2 + S_{k}'[\rho] - iik n \cdot S_{k} J + 
%  \alpha (\nabla \cdot S_{k}[J] - ik S_{k}[\rho])
%
%  This integral representation has only one scalar parameter given by
%  \alpha
%
%
%
%
%
%
%

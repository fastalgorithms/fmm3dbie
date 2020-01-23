The files are organized into the following folders

fmm_wrappers - contains wrappers which are modified
  version of the fmm wrappers and call one of the guru
  fmm interfaces under the hood

lap_wrappers - contains near quadrature extraction
  routines, and layer potential application routines
  for different integral representation for the 
  Laplace equation

helm_wrappers - contains near quadrature extraction
  routines, and layer potential application routines
  for different integral representation for the 
  Helmholtz equation

maxwell_wrappers - contains near quadrature extraction
  routines, and layer potential application routines
  for different integral representation for the 
  Maxwell equations

kernels - contains the generic kernels used 
 for determining far-field oversampling parameters
 (these are the Helmholtz single layer, double layer,
 and quadruple layer potentials)


quadratures - contains guru interfaces for evaluating
   near field quadrature corrections for compact,
   prinicpal value or hypersingular kernels for
   different quadrature schemes, estimating 
   far field oversampling for different kernels,
   estimating the near-field, and working with
   the near field structure

surface_routs - contains various routines for handling 
   functions defined on surfaces

tria_routs - routines for handling orthogonal polynomials
  on triangles and computing integrals on the standard
  triangle


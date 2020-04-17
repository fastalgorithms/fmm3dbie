Tutorial for building new wrappers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On this page, we provide a tutorial for constructing a layer potential
evaluator for 

.. math::
   
   A[\sigma] = \int_{\Gamma} G(x,y) \sigma(y) da(y)

For concreteness, we will be referring to the routines used for creating
the Helmholtz Dirichlet wrapper.


*Prerequisites* for constructing this wrapper include an FMM for
evaluating sums of the form

.. math::

   u_{i} = \sum_{j=1}^{N} G(x_{i},y_{j}) w_{j} \, , i=1,2,\ldots M \, ,

where $w_{j}$, $j=1,2,\ldots N$ are given.

To generate such a wrapper one needs to write the following minimal 
set of routines

1. :ref:`build-wrap-ker-eval`
2. :ref:`build-wrap-near`
3. :ref:`build-wrap-gurulp` (optional :ref:`build-wrap-simplelp`)
4. :ref:`build-wrap-iter`


.. _build-wrap-ker-eval:

Kernel evaluator
-----------------
A kernel evaluator for evaluating $G(x,y)$. It should have the calling
sequence.

 .. code:: fortran

   subroutine fker(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
     
Input arguments:

  - src: double precision(12)
      x,$\partial_{u} x$, $\partial_{v} x$, n
  - ndt: integer
      leading dimension of target vector
  - targ: double precision(ndt)
      target info, the first three components must be the cartesian
      coordinates of the target
  - ndd: integer
      number of double precision parameters
  - dpars: double precision(ndd)
      double precision parameters
  - ndi: integer
      number of integer parameters
  - ipars: integer(ndi)
      integer parameters
  - ndz: integer
      number of double complex parameters
  - zpars: double complex(ndz)
      complex parameters

Output parameters:
     
  - val: double complex
      value of $G(x,y)$

.. note::
  In the event the kernel is vector valued or tensor valued, we
  recommend writing a separate routine for each component for efficiency
  purposes. That being said, there are vector versions of all our base 
  routines available and the calling sequence of the vector kernel
  changes to

  .. code:: fortran

    subroutine fker_vec(nd,src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)

  The only difference being now the output $val$ is an array of size nd.

Example::

  src/kernels/helm_kernels.f90
  subroutine h3d_slp

.. _build-wrap-near:

Near quadrature correction wrapper
-----------------------------------
   
Should be analogous to the `Helmholtz Dirichlet quadrature generation
wrapper <helm_wrappers.html#getnearquad-helm-comb-dir>`__. 

Subroutine summary
===================

This routine internally calls a guru version wrapper for generating
near quadrature corrections depending on whether the kernel is
compact, pv, or hypersingular::

  call getnearquad_ggq_compact_guru(..)
  call getnearquad_ggq_pv_guru(..)


.. _build-wrap-gurulp:

Guru layer potential evaluator routine
---------------------------------------

Should be analogous to the `Guru interface for helmholtz combined field
dirichlet layer potential evaluator
<helm_wrappers.html#lpcomp-helm-comb-dir-addsub>`_.
The subroutine given near field information, precomputed quadrature correction, 
oversampled geometry information for far-field, and a density, evaluates
the layer potential.


Example::
  
  src/helm_wrappers/helm_comb_dir.f
  subroutine lpcomp_helm_comb_dir_addsub(..)

Subroutine summary
===================

1. Oversample the density::
   
     call oversample_fun_surf(..)

2. Call the fmm::
 
     call hfmm3d(...)

3. Add in precomputed quadrature correction

4. Subtract corresponding contributions computed via the FMM

.. note::
   
   Care must be taken in using the appropriate source grid in the add
   and subtract step of evaluating the layer potentila. 
   The FMM computes interactions from the oversampled grid, while the
   near quadrature correction acts on the original discretization.  

.. _build-wrap-simplelp:

Simple layer potential evaluator routine
-----------------------------------------

Should be analogous to the `Helmholtz combined field
dirichlet layer potential evaluator
<helm_wrappers.html#lpcomp-helm-comb-dir>`_.
The subroutine given a density, and a collection of targets, evaluates
the layer potential.


Example::
  
  src/helm_wrappers/helm_comb_dir.f
  subroutine lpcomp_helm_comb_dir(..)

Subroutine summary
===================

1. Determine near field::

     call get_centroid_rads(..): Compute centroid and bounding sphere radii
     call get_rfacs(..): Estimates near field parameters
     call findnearmem(..): Estimates memory requirements for row sparse compressed structure
     call findnear(..): Computes the row sparse compressed structure 
     call get_iquad_rsc(..): estimate iquad pointer 

2. Oversample the geometry::
 
     call get_far_order(..): estimate required upsampling
     call oversample_geom(..): oversample the surface information
     call get_qwts(..): get oversampledquadrature weights 
     
3. Compute near quadrature correction::    

     call getnearquad_helm_comb_dir(..) 

4. Evaluate layer potential using guru interface::

     call lpcomp_helm_comb_dir_addsub(..) 

.. _build-wrap-iter:

Iterative solver routine
------------------------

Should be analogous to the `Helmholtz combined field
dirichlet problem iterative solver
<helm_wrappers.html#helm-comb-dir-iter-solver>`_.
The subroutine given a density, and a collection of targets, evaluates
the layer potential.


Example::
  
  src/helm_wrappers/helm_comb_dir.f
  subroutine lpcomp_helm_comb_dir(..)

Subroutine summary
===================

1. Setup targets::
     
     call patch_id_uvs(..)

2. Determine near field::
   
     call get_centroid_rads(..): Compute centroid and bounding sphere radii
     call get_rfacs(..): Estimates near field parameters
     call findnearmem(..): Estimates memory requirements for row sparse compressed structure
     call findnear(..): Computes the row sparse compressed structure 
     call get_iquad_rsc(..): estimate iquad pointer 

3. Oversample the geometry::
 
     call get_far_order(..): estimate required upsampling
     call oversample_geom(..): oversample the surface information
     call get_qwts(..): get oversampledquadrature weights 
     
4. Compute near quadrature correction::    

     call getnearquad_helm_comb_dir(..) 

5. Compute solution using GMRES (avoiding low-threshold stagnation)  


.. note::

   In the GMRES solver for iterating $(\alpha I + K)x=y$, you construct the
   Krylov subspace for $K^{n}y$, where $b$ is the given data, and
   manually add in the identity term. This is handled by the variable
   $\alpha=$ ``zid``.

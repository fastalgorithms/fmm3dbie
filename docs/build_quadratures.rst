Tutorial for building new quadrature schemes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On this page, we provide a tutorial for implementing a new quadrature
approach.

For concreteness, we will be referring to the Generalized Gaussian
quadrautre + adaptive integration approach for compact kernels::

  src/quadratures/ggq-quads.f
  subroutine getnearquad_ggq_compact_guru

Subroutine summary
-------------------

1. Convert row sparse compressed format to column sparse compressed
   format (This is needed in case you wish to amortize precomputed
   source information for different targets close to the same patch)::
     
     call rsc_to_csc(..)

2. Compute centroids and bounding sphere radii (this is needed in case
   you wish to use oversampled smooth quadrature for targets that are
   not too close to the patch)::

     call get_centroid_rads(..)

3. Get adaptive integration parameters (these have been optimized for
   adaptive integration for the standard kernels)::

     call get_quadparams_adap

4. Compute the integrals:
   
   - split the near field into list of
     targets to be handled via adaptive integration, and list of targets
     to be handled via oversampled quadrature

   .. code:: fortran

     call get_near_far_split_pt(..)
   
   - For a given patch, comptute integrals of all
     orthogonal polynomials as density on the given patch at all near-near
     targets using adaptive integration (later on we
     use the vals-to-coefs matrix to get the near quadrature
     correction)

   .. code:: fortran    

     call ctriaints(..) 

   - Oversampled quadrature for targets not too close (near-far targets)

   .. code:: fortran
   
     call ctriaints_wnodes(..): 

   - Specialized quadratures for self interaction

   .. code:: fortran

     call get_ggq_self_quad_pt(..)

.. note::
   A remark for quadrature approaches like QBX, where the kernel itself
   depends on the target location for example:
   
   .. math::
      
      G(t,s) = \sum_{\ell=0}^{p} f_{\ell}(|t-c_{t}|) g_{\ell}(|s-c_{t}|) \, , 
      
   One could send in the parameters,  :math:`c_{t},f_{\ell}(|t-c_{t}|)` as a
   collection of additional target parameters in the ``targ`` array.  

.. note::
   To contribute a quadrature method to the library, update at least one
   of the wrappers in ``src/helm_wrappers``, or ``src/lap_wrappers``, or
   ``src/maxw_wrappers`` or ``src/stok_wrappers`` (preferably all of
   them), and submit a pull request

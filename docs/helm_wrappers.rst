.. _hwrap:

Helmholtz solvers
==================

- :ref:`helm-dir`

.. _helm-dir:

Dirichlet boundary value problem
**********************************

For solving the Helmholtz double layer potential, we use the combined
field representation:

.. math::

   u = \alpha \mathcal{S}_{k}[\sigma] + \beta \mathcal{D}_{k}[\sigma] \,.

We have the following user callable routines:

- :ref:`lpcomp-helm-comb-dir`
- :ref:`helm-comb-dir-iter-solver`
- :ref:`helm-dir-adv`

.. _lpcomp-helm-comb-dir:

Layer potential evaluator
~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluate dirichlet data for Helmholtz equation using combined field
integral equation

.. code:: fortran

   subroutine lpcomp_helm_comb_dir(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,ndtarg,targs,ipatch_id,uvs_targ,eps,zpars,sigma,pot)

This subroutine evaluates the layer potential for the representation

.. math ::

   u = 4 \pi (\alpha \mathcal{S}_{k}[\sigma] + \beta D_{k}[\sigma] )

Note the additional scaling of $4\pi$ to be consistent with the FMM3D
libraries. For targets on the boundary, this routine only computes the
principal value part, the identity term corresponding to the jump in the
double layer is not included in the layer potential.

Input arguments:

    - npatches: integer
        number of patches
    - norders: integer(npatches)
        order of discretization on each patch
    - ixyzs: integer(npatches+1)
        ixyzs(i) denotes the starting location in srccoefs, and
        srcvals array where information for patch i begins
    - iptype: integer(npatches)
        type of patch
    - npts: integer
        total number of points on the boundary
    - srccoefs: double precision (9,npts)
        koornwinder exapansion coefficients for x, $\partial_{u} x$, and
        $\partial_{v} x$
    - srcvals: double precision (12,npts)
        x,$\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
        discretization nodes
    - ndtarg: integer
        leading dimension of target array
    - ntarg: integer
        number of targets
    - targs: double precision (ndtarg,ntarg)
        target information. Note the first three entries must be the
        cartesian components of the target
    - ipatch_id: integer(ntarg)
         id of patch of target i, id= -1/0 if target is off-surface
    - uvs_targ: double precision (2,ntarg)
         local uv coordinates on patch, if target is on surface
    - eps: double precision
         precision requested
    - zpars: double complex(3)
         kernel parameters, zpars(1)=k, zpars(2)=$\alpha$,
         zpars(3)=$\beta$
    - sigma: double complex(npts)
         density for layer potential
    
Output arguments:
    - pot: double complex(ntarg)
         layer potential evaluated at target points

.. container:: rttext

   `Back to Dirichlet solvers <helm_wrappers.html#helm-dir>`__

.. container:: rttext

   `Back to top <helm_wrappers.html#hwrap>`__

.. _helm-comb-dir-iter-solver:

Iterative solver
~~~~~~~~~~~~~~~~~~~~~~~~~~


Evaluate dirichlet data for Helmholtz equation using combined field
integral equation

.. code:: fortran

   subroutine helm_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres,niter,errs,
   rres,sigma)

This subroutine solves the interior or exterior Dirichlet boundary value
problem using the combined field integral representation

.. math ::

   u = 4 \pi (\alpha \mathcal{S}_{k}[\sigma] + \beta D_{k}[\sigma] )

The integral equation on the boundary is given by

.. math::
   
   f = \pm 2\pi \beta \sigma + 4\pi (\alpha \mathcal{S}_{k}[\sigma] +
   \beta D_{k}[\sigma]) \, ,

where f is the boundary data, and sign of the constant is positive for
the exterior problem, and negative for the interior problem. 

Note the additional scaling of $4\pi$ to be consistent with the FMM3D
libraries. 

Input arguments:

    - npatches: integer
        number of patches
    - norders: integer(npatches)
        order of discretization on each patch
    - ixyzs: integer(npatches+1)
        ixyzs(i) denotes the starting location in srccoefs, and
        srcvals array where information for patch i begins
    - iptype: integer(npatches)
        type of patch
    - npts: integer
        total number of points on the boundary
    - srccoefs: double precision (9,npts)
        koornwinder exapansion coefficients for x, $\partial_{u} x$, and
        $\partial_{v} x$
    - srcvals: double precision (12,npts)
        x,$\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
        discretization nodes
    - eps: double precision
         precision requested for computing quadrature and fmm tolerance
    - zpars: double complex(3)
         kernel parameters, zpars(1)=k, zpars(2)=$\alpha$,
         zpars(3)=$\beta$
    - ifinout: integer
         flag fpr interior or exterior problems (normals assumed to be
         pointing in the exertior of the region)
         -  ifinout = 0, interior problem
         -  ifinout = 1, exterior problem
    - rhs: double complex(npts)
         boundary data f
    - eps_gmres: double precision
         gmres tolerance requested
    - numit: integer
         max number of gmres iterations
    
Output arguments:
    - niter: integer
         number of gmres iterations required for relative residual
    - errs: double precision(1:niter)
         relative residual as a function of iteration number
    - rres: double precision
         relative residual of the computed solution
    - sigma: double complex(npts)
         solution of integral equation
    

.. container:: rttext

   `Back to Dirichlet solvers <helm_wrappers.html#helm-dir>`__

.. container:: rttext

   `Back to top <helm_wrappers.html#hwrap>`__

.. _helm-dir-adv:

Advanced wrappers
~~~~~~~~~~~~~~~~~~

We also the following advanced user interfaces

- :ref:`getnearquad-helm-comb-dir`
- :ref:`lpcomp-helm-comb-dir-addsub`
- :ref:`helm-fds-routs`

.. _getnearquad-helm-comb-dir:

Near quadrature correction generation
---------------------------------------


Evaluate near quadrature correction and store in row sparse
compressed format 

.. code:: fortran

   subroutine getnearquad_helm_comb_dir(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps,zpars,
   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)

This subroutine generates the near field quadrature for the
representation

.. math ::

   u = 4 \pi (\alpha \mathcal{S}_{k}[\sigma] + \beta D_{k}[\sigma] ) \,,

where the near field is specified by the user in row sparse compressed
format

Note the additional scaling of $4\pi$ to be consistent with the FMM3D
libraries. 

Input arguments:

    - npatches: integer
        number of patches
    - norders: integer(npatches)
        order of discretization on each patch
    - ixyzs: integer(npatches+1)
        ixyzs(i) denotes the starting location in srccoefs, and
        srcvals array where information for patch i begins
    - iptype: integer(npatches)
        type of patch
    - npts: integer
        total number of points on the boundary
    - srccoefs: double precision (9,npts)
        koornwinder exapansion coefficients for x, $\partial_{u} x$, and
        $\partial_{v} x$
    - srcvals: double precision (12,npts)
        x,$\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
        discretization nodes
    - ndtarg: integer
        leading dimension of target array
    - ntarg: integer
        number of targets
    - targs: double precision (ndtarg,ntarg)
        target information. Note the first three entries must be the
        cartesian components of the target
    - ipatch_id: integer(ntarg)
         id of patch of target i, id= -1/0 if target is off-surface
    - uvs_targ: double precision (2,ntarg)
         local uv coordinates on patch, if target is on surface
    - eps: double precision
         precision requested for computing quadrature and fmm tolerance
    - zpars: double complex(3)
         kernel parameters, zpars(1)=k, zpars(2)=$\alpha$,
         zpars(3)=$\beta$
    - iquadtype: integer
         quadrature type

         - iquadtype = 1: for generalized gaussian quadrature+adaptive
           integration
    - nnz: integer
         number of non-zero target-patch interactions
    - row_ptr: integer(ntarg+1) 
         row_ptr(i) is the starting location in the
         col-ind list where the list of patches in the near field of target
         i start. If row_ptr(i) <= j < row_ptr(i+1), then target(i) and
         patch col_ind(j) are in the near-field of each other
    - col_ind: integer(nnz)
         list of patches which interact with the targets. See `example
         <near_reps.html#nearreps-exmp>`_
    - iquad: integer(nnz)
          iquad(i) is the location in the quadrature array where the
          matrix entries corresponding to the interaction of target i,
          and patch col_ind(j) start in wnear array
    - rfac0: double precision
          radius parameter for switiching between adaptive integration
          and oversampled smooth quadrature in quadrature generation
          routine (note this is different from the radius parameter 
          used for identifying the near field)
    - nquad: integer
          size of wnear array

Output arguments:
    - wnear: double complex(nquad)
          near field quadrature correction


.. container:: rttext

   `Back to top <helm_wrappers.html#hwrap>`__

.. _lpcomp-helm-comb-dir-addsub:

Guru interface for layer potential evaluation
----------------------------------------------

Guru interface for evaluating dirichlet data for Helmholtz equation 
using combined field integral equation

.. code:: fortran

   subroutine lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,ndtarg,targs,ipatch_id,uvs_targ,eps,zpars,iquadtype,
   nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,nptso,ixyzso,
   srcover,wover,pot)

This subroutine evaluates the layer potential for the representation

.. math ::

   u = 4 \pi (\alpha \mathcal{S}_{k}[\sigma] + \beta D_{k}[\sigma] )

Note the additional scaling of $4\pi$ to be consistent with the FMM3D
libraries. For targets on the boundary, this routine only computes the
principal value part, the identity term corresponding to the jump in the
double layer is not included in the layer potential.

Input arguments:

    - npatches: integer
        number of patches
    - norders: integer(npatches)
        order of discretization on each patch
    - ixyzs: integer(npatches+1)
        ixyzs(i) denotes the starting location in srccoefs, and
        srcvals array where information for patch i begins
    - iptype: integer(npatches)
        type of patch
    - npts: integer
        total number of points on the boundary
    - srccoefs: double precision (9,npts)
        koornwinder exapansion coefficients for x, $\partial_{u} x$, and
        $\partial_{v} x$
    - srcvals: double precision (12,npts)
        x,$\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
        discretization nodes
    - ndtarg: integer
        leading dimension of target array
    - ntarg: integer
        number of targets
    - targs: double precision (ndtarg,ntarg)
        target information. Note the first three entries must be the
        cartesian components of the target
    - ipatch_id: integer(ntarg)
         id of patch of target i, id= -1/0 if target is off-surface
    - uvs_targ: double precision (2,ntarg)
         local uv coordinates on patch, if target is on surface
    - eps: double precision
         precision requested
    - zpars: double complex(3)
         kernel parameters, zpars(1)=k, zpars(2)=$\alpha$,
         zpars(3)=$\beta$
    - iquadtype: integer
         quadrature type

         - iquadtype = 1: for generalized gaussian quadrature+adaptive
           integration
    - nnz: integer
         number of non-zero target-patch interactions
    - row_ptr: integer(ntarg+1) 
         row_ptr(i) is the starting location in the
         col-ind list where the list of patches in the near field of target
         i start. If row_ptr(i) <= j < row_ptr(i+1), then target(i) and
         patch col_ind(j) are in the near-field of each other
    - col_ind: integer(nnz)
         list of patches which interact with the targets. See `example
         <near_reps.html#nearreps-exmp>`_
    - iquad: integer(nnz)
          iquad(i) is the location in the quadrature array where the
          matrix entries corresponding to the interaction of target i,
          and patch col_ind(j) start in wnear array
    - nquad: integer
          size of wnear array
    - wnear: double complex(nquad)
          near field quadrature correction
    - sigma: double complex(npts)
         density for layer potential
    - novers: integer(npatches)
        order of oversampled discretization on each patch
    - nptso: integer
        number of oversampled discretization nodes
    - ixyzso: integer(npatches+1)
        ixyzso(i) denotes the starting location in srcover, and
        wover array where information for patch i begins
    - srcover: double precision (12,nptso)
        x,$\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
        discretization nodes
    - wover: double precision(nptso)
        oversampled quadrature weights
    
Output arguments:
    - pot: double complex(ntarg)
         layer potential evaluated at target points

.. container:: rttext

   `Back to top <helm_wrappers.html#hwrap>`__

.. _helm-fds-routs:

Fast direct solver routines
----------------------------

Under construction

.. container:: rttext

   `Back to top <helm_wrappers.html#hwrap>`__


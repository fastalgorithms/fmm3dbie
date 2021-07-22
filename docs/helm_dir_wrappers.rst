.. _hdirwrap:

Helmholtz dirichlet wrappers
=============================

- `Dirichlet iterative solver <helm_dir_wrappers.html#helm-comb-dir-iter-solver>`__
- `Dirichlet post-processor <helm_dir_wrappers.html#lpcomp-helm-comb-dir>`__
- `Advanced dirichlet wrappers <helm_dir_wrappers.html#helm-dir-adv>`__ 

.. _helm-comb-dir-iter-solver:

Iterative solver
~~~~~~~~~~~~~~~~~~~~~~~~~~

Solve the dirichlet boundary value problem for the Helmholtz equation 
using combined field integral equation

.. code:: fortran

   subroutine helm_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres,niter,errs,
   rres,sigma)

.. include:: raws/helm_dir/helm-comb-dir-solver.raw

.. container:: rttext

   `Back to Helmholtz solvers <helm_wrappers.html#hwrap>`__

.. container:: rttext

   `Back to top <helm_dir_wrappers.html#hdirwrap>`__


.. _lpcomp-helm-comb-dir:

Dirichlet problem post-processor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluate potential for the Helmholtz equation using combined field
integral equation

.. code:: fortran

   subroutine lpcomp_helm_comb_dir(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,ndtarg,targs,ipatch_id,uvs_targ,eps,zpars,sigma,pot)

.. include:: raws/helm_dir/lpcomp-helm-comb-dir.raw   

.. container:: rttext

   `Back to Helmholtz solvers <helm_wrappers.html#hwrap>`__

.. container:: rttext

   `Back to top <helm_dir_wrappers.html#hdirwrap>`__

.. _helm-dir-adv:

Advanced wrappers
~~~~~~~~~~~~~~~~~~

We also the following advanced user interfaces

- :ref:`getnearquad-helm-comb-dir`
- :ref:`lpcomp-helm-comb-dir-addsub`

.. _getnearquad-helm-comb-dir:

Near quadrature correction generation
---------------------------------------


Evaluate near quadrature correction and store in row sparse
compressed format 

.. code:: fortran

   subroutine getnearquad_helm_comb_dir(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps,zpars,
   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)

.. include:: raws/helm_dir/getnearquad-helm-comb-dir.raw


.. container:: rttext

   `Back to Helmholtz solvers <helm_wrappers.html#hwrap>`__

.. container:: rttext

   `Back to top <helm_dir_wrappers.html#hdirwrap>`__

.. _lpcomp-helm-comb-dir-addsub:

Guru interface for layer potential evaluation
----------------------------------------------

Guru interface for evaluating the potential for the Helmholtz equation 
using combined field integral equation

.. code:: fortran

   subroutine lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,ndtarg,targs,ipatch_id,uvs_targ,eps,zpars,
   nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,nptso,ixyzso,
   srcover,wover,pot)

.. include:: raws/helm_dir/lpcomp-helm-comb-dir-addsub.raw

.. container:: rttext

   `Back to Helmholtz solvers <helm_wrappers.html#hwrap>`__

.. container:: rttext

   `Back to top <helm_dir_wrappers.html#hdirwrap>`__



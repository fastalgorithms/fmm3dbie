.. _hneuwrap:

Helmholtz Neumann wrappers
=============================

- `Neumann iterative solver <helm_neu_wrappers.html#helm-rpcomb-neu-iter-solver>`__
- `Neumann post-processor <helm_neu_wrappers.html#lpcomp-helm-rpcomb-dir>`__

.. _helm-rpcomb-neu-iter-solver:

Iterative solver
~~~~~~~~~~~~~~~~~~~~~~~~~~

Solve neumann boundary value problem for the Helmholtz equation 
using the right preconditioned combined field integral equation

.. code:: fortran

   subroutine helm_rpcomb_neu_solver(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres,niter,errs,
   rres,sigma,siksigma)

.. include:: raws/helm_neu/helm-rpcomb-neu-solver.raw

.. container:: rttext

   `Back to Helmholtz solvers <helm_wrappers.html#hwrap>`__

.. container:: rttext

   `Back to top <helm_neu_wrappers.html#hneuwrap>`__


.. _lpcomp-helm-rpcomb-dir:

Neumann problem post-processor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluate potential for the Helmholtz equation using right preconditioned 
combined field integral equation

.. code:: fortran

   subroutine lpcomp_helm_rpcomb_dir(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,ndtarg,targs,ipatch_id,uvs_targ,eps,zpars,sigma,
   siksigma,pot)

.. include:: raws/helm_neu/lpcomp-helm-rpcomb-dir.raw   

.. container:: rttext

   `Back to Helmholtz solvers <helm_wrappers.html#hwrap>`__

.. container:: rttext

   `Back to top <helm_neu_wrappers.html#hneuwrap>`__



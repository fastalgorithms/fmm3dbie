.. _himpwrap:

Helmholtz Impedance wrappers
=============================


- `Impedance iterative solver <helm_imp_wrappers.html#helm-rpcomb-imp-iter-solver>`__
- `Impedance post-processor <helm_imp_wrappers.html#lpcomp-helm-rpcomb-dir-imp>`__

.. _helm-rpcomb-imp-iter-solver:

Iterative solver
~~~~~~~~~~~~~~~~~~~~~~~~~~

Solve impedance boundary value problem for the Helmholtz equation 
using the right preconditioned combined field integral equation

.. code:: fortran

   subroutine helm_rpcomb_imp_solver(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,eps,zpars,zlams,numit,ifinout,rhs,eps_gmres,niter,errs,
   rres,sigma,siksigma)

.. include:: raws/helm_imp/helm-rpcomb-imp-solver.raw

.. container:: rttext

   `Back to Helmholtz solvers <helm_wrappers.html#hwrap>`__

.. container:: rttext

   `Back to top <helm_imp_wrappers.html#himpwrap>`__


.. _lpcomp-helm-rpcomb-dir-imp:

Impedance problem post-processor
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

   `Back to top <helm_imp_wrappers.html#himpwrap>`__



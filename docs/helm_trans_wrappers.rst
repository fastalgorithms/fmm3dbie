.. _htranswrap:

Helmholtz Transmission wrappers
=================================

- `Transmission iterative solver <helm_trans_wrappers.html#helm-comb-trans-iter-solver>`__
- `Transmission post-processor <helm_trans_wrappers.html#lpcomp-comb-split-dir>`__

.. _helm-comb-trans-iter-solver:

Iterative solver
~~~~~~~~~~~~~~~~~~~~~~~~~~

Solve impedance boundary value problem for the Helmholtz equation 
using the right preconditioned combined field integral equation

.. code:: fortran

   subroutine helm_comb_trans_solver(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,eps,zpars,numit,rhs,eps_gmres,niter,errs,
   rres,sigma,siksigma)

.. include:: raws/helm_trans/helm-comb-trans-solver.raw

.. container:: rttext

   `Back to Helmholtz solvers <helm_wrappers.html#hwrap>`__

.. container:: rttext

   `Back to top <helm_trans_wrappers.html#htranswrap>`__


.. _lpcomp-helm-comb-split-dir:

Transmission problem post-processor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluate potential for the Helmholtz equation using  
combined field integral equation where a different density is given
for the single and double layer potentials

.. code:: fortran

   subroutine lpcomp_helm_comb_split_dir(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,ndtarg,targs,ipatch_id,uvs_targ,eps,zpars,sigma,
   pot)

.. include:: raws/helm_trans/lpcomp-helm-comb-split-dir.raw   

.. container:: rttext

   `Back to Helmholtz solvers <helm_wrappers.html#hwrap>`__

.. container:: rttext

   `Back to top <helm_trans_wrappers.html#htranswrap>`__



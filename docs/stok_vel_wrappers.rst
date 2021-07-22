.. _svelwrap:

Stokes velocity wrappers
=============================

- `Dirichlet iterative solver <stok_vel_wrappers.html#stok-comb-vel-iter-solver>`__
- `Dirichlet post-processor <stok_vel_wrappers.html#lpcomp-stok-comb-vel>`__

.. _stok-comb-vel-iter-solver:

Iterative solver
~~~~~~~~~~~~~~~~~~~~~~~~~~

Solve the velocity boundary value problem for the Stokes equation 
using combined field integral equation

.. code:: fortran

   subroutine stok_comb_vel_solver(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,eps,dpars,numit,ifinout,rhs,eps_gmres,niter,errs,
   rres,sigma)

.. include:: raws/stok_vel/stok-comb-vel-solver.raw

.. container:: rttext

   `Back to Stokes solvers <stok_wrappers.html#swrap>`__

.. container:: rttext

   `Back to top <stok_vel_wrappers.html#svelwrap>`__


.. _lpcomp-stok-comb-vel:

Velocity problem post-processor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluate velocity for the Stokes equation using combined field
integral equation

.. code:: fortran

   subroutine lpcomp_stok_comb_vel(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,ndtarg,targs,ipatch_id,uvs_targ,eps,dpars,sigma,pot)

.. include:: raws/stok_vel/lpcomp-stok-comb-vel.raw   

.. container:: rttext

   `Back to Stokes solvers <stok_wrappers.html#swrap>`__

.. container:: rttext

   `Back to top <stok_vel_wrappers.html#svelwrap>`__

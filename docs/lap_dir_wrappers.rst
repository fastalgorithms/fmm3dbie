.. _ldirwrap:

Laplace dirichlet wrappers
=============================

- `Dirichlet iterative solver <lap_dir_wrappers.html#lap-comb-dir-iter-solver>`__
- `Dirichlet post-processor <lap_dir_wrappers.html#lpcomp-lap-comb-dir>`__

.. _lap-comb-dir-iter-solver:

Iterative solver
~~~~~~~~~~~~~~~~~~~~~~~~~~

Solve the dirichlet boundary value problem for the Laplace equation 
using combined field integral equation

.. code:: fortran

   subroutine lap_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,eps,dpars,numit,ifinout,rhs,eps_gmres,niter,errs,
   rres,sigma)

.. include:: raws/lap_dir/lap-comb-dir-solver.raw

.. container:: rttext

   `Back to Laplace solvers <lap_wrappers.html#lwrap>`__

.. container:: rttext

   `Back to top <lap_dir_wrappers.html#ldirwrap>`__


.. _lpcomp-lap-comb-dir:

Dirichlet problem post-processor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluate potential for the Laplace equation using combined field
integral equation

.. code:: fortran

   subroutine lpcomp_lap_comb_dir(npatches,norders,ixyzs,iptype,npts,
   srccoefs,srcvals,ndtarg,targs,ipatch_id,uvs_targ,eps,dpars,sigma,pot)

.. include:: raws/lap_dir/lpcomp-lap-comb-dir.raw   

.. container:: rttext

   `Back to Laplace solvers <lap_wrappers.html#lwrap>`__

.. container:: rttext

   `Back to top <lap_dir_wrappers.html#ldirwrap>`__

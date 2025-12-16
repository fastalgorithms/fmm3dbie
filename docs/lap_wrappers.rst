.. _lwrap:

Laplace solvers
==================

- :ref:`lap-dir`

.. _lap-dir:

Dirichlet boundary value problem
*******************************************

The dirichlet problem for the Laplace equation is given by

.. math::

   \Delta u &= 0 \quad \mbox{ in } \Omega \, , \\
   u &= f \quad \mbox{ on } \Gamma \, .

   
For solving the above problem, we represent $u$ in $\Omega$ using
the combined field representation:

.. math::

   u = \alpha \mathcal{S}_{0}[\sigma] + \beta \mathcal{D}_{0}[\sigma] \,.

We have the following user callable routines:

- `Dirichlet iterative solver <lap_dir_wrappers.html#lap-comb-dir-iter-solver>`__
- `Dirichlet post-processor <lap_dir_wrappers.html#lpcomp-lap-comb-dir>`__

To see a demo of the use of these interfaces see
``examples/laplace/lap_dir_iter_example.f``. 
This script can be run using ``make -f lap_dir_iter_example.make`` in the
``examples/laplace`` folder.

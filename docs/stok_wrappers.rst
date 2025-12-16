.. _swrap:

Stokes solvers
==================

- :ref:`stok-vel`

.. _stok-vel:

Velocity boundary value problem
*******************************************

The velocity problem for the Stokes equation is given by

.. math::

   \Delta u &= \nabla p \quad \mbox{ in } \Omega \, , \\
   \nabla \cdot  u &= 0 \quad \mbox{ in } \Omega \, , \\
   u &= f \quad \mbox{ on } \Gamma \, .

   
For solving the above problem, we represent $u$ in $\Omega$ using
the combined field representation:

.. math::

   u = \alpha \mathcal{S}^{\textrm{stok}}[\sigma] + \beta \mathcal{D}^{\textrm{stok}}[\sigma] \,.

We have the following user callable routines:

- `Velocity iterative solver <stok_vel_wrappers.html#stok-comb-vel-iter-solver>`__
- `Velocity post-processor <stok_vel_wrappers.html#lpcomp-stok-comb-vel>`__

To see a demo of the use of these interfaces see
``examples/stokes/stok_vel_iter_example.f``. 
This script can be run using ``make -f stok_vel_iter_example.make`` in the
``examples/stokes`` folder.

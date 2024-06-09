.. _hwrap:

Helmholtz solvers
==================

- :ref:`helm-dir`
- :ref:`helm-neu`
- :ref:`helm-imp`
- :ref:`helm-trans`

.. _helm-dir:

Dirichlet boundary value problem
*******************************************

The dirichlet problem for the Helmholtz equation is given by

.. math::

   (\Delta + k^2) u &= 0 \quad \mbox{ in } \Omega \, , \\
   u &= f \quad \mbox{ on } \Gamma \, .

   
For solving the above problem, we represent $u$ in $\Omega$ using
the combined field representation:

.. math::

   u = \alpha \mathcal{S}_{k}[\sigma] + \beta \mathcal{D}_{k}[\sigma] \,.

We have the following user callable routines:

- `Dirichlet iterative solver <helm_dir_wrappers.html#helm-comb-dir-iter-solver>`__
- `Dirichlet post-processor <helm_dir_wrappers.html#lpcomp-helm-comb-dir>`__

To see a demo of the use of these interfaces see
``examples/helmholtz/helm_dir_iter_example.f``. 
This script can be run using ``make -f helm_dir_iter_example.make`` in the
``examples/helmholtz`` folder.

For advanced Dirichlet wrappers see:

- `Advanced dirichlet wrappers <helm_dir_wrappers.html#helm-dir-adv>`__ 

.. _helm-neu:

Neumann boundary value problem
*******************************************

The Neumann problem for the Helmholtz equation is given by

.. math::

   (\Delta + k^2) u &= 0 \quad \mbox{ in } \Omega \, , \\
   \frac{\partial u}{\partial n} &= f \quad \mbox{ on } \Gamma \, .

   
For solving the above problem, we represent $u$ in $\Omega$ using
the right preconditioned combined field representation:

.. math::

   u = \mathcal{S}_{k}[\sigma] + i\alpha \mathcal{D}_{k}[S_{i|k|}[\sigma]] \,.

We have the following user callable routines:

- `Neumann iterative solver <helm_neu_wrappers.html#helm-rpcomb-neu-iter-solver>`__
- `Neumann post-processor <helm_neu_wrappers.html#lpcomp-helm-rpcomb-dir>`__

To see a demo of the use of these interfaces see
``examples/helmholtz/helm_neu_iter_example.f``. 
This script can be run using ``make -f helm_neu_iter_example.make`` in the
``examples/helmholtz`` folder.

.. _helm-imp:

Impedance boundary value problem
*******************************************

The Impedance problem for the Helmholtz equation is given by

.. math::

   (\Delta + k^2) u &= 0 \quad \mbox{ in } \Omega \, , \\
   \frac{\partial u}{\partial n} + ik \lambda u &= f \quad \mbox{ on } \Gamma \, .

   
For solving the above problem, we represent $u$ in $\Omega$ using
the right preconditioned combined field representation:

.. math::

   u = \mathcal{S}_{k}[\sigma] + i\alpha \mathcal{D}_{k}[S_{i|k|}[\sigma]] \,.

We have the following user callable routines:

- `Impedance iterative solver <helm_imp_wrappers.html#helm-rpcomb-imp-iter-solver>`__
- `Impedance post-processor <helm_imp_wrappers.html#lpcomp-helm-rpcomb-dir-imp>`__

To see a demo of the use of these interfaces see
``examples/helmholtz/helm_imp_iter_example.f``. 
This script can be run using ``make -f helm_imp_iter_example.make`` in the
``examples/helmholtz`` folder.


.. _helm-trans:

Transmission boundary value problem
*******************************************

The Impedance problem for the Helmholtz equation is given by

.. math::

   (\Delta + k_{1}^2) u_{1} &= 0 \quad \mbox{ in } \Omega \, , \\
   (\Delta + k_{0}^2) u_{0} &= 0 \quad \mbox{ in } \Omega^c \, , \\
   \alpha_{0} u_{0} - \alpha_{1} u_{1} &= f \quad \mbox{ on } \Gamma \, , \\
   \beta_{0} \frac{\partial u_{0}}{\partial n} - \beta_{1} \frac{\partial u_{1}}{\partial n} &= g \quad \mbox{ on } \Gamma \, .

where $k_{0, k_{1}$, are the interior and exterior wave numbers respectively
   
For solving the above problem, we represent $u$ using the combined field
representation:

.. math::

   u_{1} &= \frac{1}{\beta_{1}} \left( ik_{1} \mathcal{S}_{k_{1}}[\lambda] + \mathcal{D}_{k_{1}}[\rho] \right) \,, \\
   u_{0} &= \frac{1}{\beta_{0}} \left( ik_{0} \mathcal{S}_{k_{0}}[\lambda] + \mathcal{D}_{k_{0}}[\rho] \right) \,.

We have the following user callable routines:

- `Transmission iterative solver <helm_trans_wrappers.html#helm-comb-trans-iter-solver>`__
- `Transmission post-processor <helm_trans_wrappers.html#lpcomp-helm-comb-split-dir>`__

To see a demo of the use of these interfaces see
``examples/helmholtz/helm_trans_iter_example.f``. 
This script can be run using ``make -f helm_trans_iter_example.make`` in the
``examples/helmholtz`` folder.

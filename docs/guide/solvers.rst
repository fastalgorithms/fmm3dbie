

.. role:: matlab(code)
   :language: matlab   

Solvers and evaluators
=======================

In MATLAB, there is a guru solver (:matlab:`solver`) routine, and an evaluator routine (:matlab:`eval_fields`).

The solvers, and evaluators can also be accessed using 
``<pde>.solver`` or ``<pde>.<bc>.solver``, and ``<pde>.eval`` or ``<pde>.<bc>.eval``
respectively.
For example, for the Helmholtz Dirichlet solver, you can use 
:matlab:`solver` or :matlab:`helm3d.solver` or :matlab:`helm3d.dirichlet.solver`.
Similarly, for the Laplace Neumann post-processor, you can use
:matlab:`eval_fields` or :matlab:`lap3d.eval` or :matlab:`lap3d.neumann.solver`.

- <pde>: partial differential equation being solved
  
  - lap3d: Laplace solvers
  - helm3d: Helmholtz solvers
  - stok3d: Stokes solvers
  - em3d: Maxwell solvers

- <bc>: boundary condition

  - dirichlet: Dirichlet boundary conditions (available with lap3d, and helm3d)
  - neumann: Neumann boundary conditions (available with lap3d, and helm3d)
  - impedance: Impedance boundary conditions (available with helm3d)
  - transmission: Transmission boundary conditions (available with helm3d)
  - velocity: Velocity boundary conditions (available with stok3d)
  - pec: Perfect electric condudctor boundary conditions (available with em3d)


solver
------

The following boundary value problems can be solved using the :matlab:`solver` function:

.. include:: ../../matlab/solver.m
   :literal:
   :code: matlab:
   :end-before: %---------------


eval_fields
------------

The following integral representations can be evaluated using the :matlab:`eval_fields` function:

.. include:: ../../matlab/eval_fields.m
   :literal:
   :code: matlab:
   :end-before: %---------------


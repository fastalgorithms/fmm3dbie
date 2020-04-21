Boundary value problem solvers
===============================

For each of the solver in this repository, we have the following 
user callable routines::

    - lpcomp_<pde>_<rep>_<bvp>: evaluate layer potential
    - <pde>_<rep>_<bvp>_solver: iterative solver

- <pde>: PDE being solved.

    - helm: Helmholtz
    - lap: Laplace
    - stok: Stokes
    - maxw: Maxwell

- <rep>: Layer potential representation for solution. Some examples include

    - s: single layer
    - d: double layer
    - comb: combined field representation
    - efie: electric field integral representation
    - mfie: magnetic field integral representation
    - dpie: decoupled potential integral representation

- <bvp>: Boundary value problem. Some examples include

    - dir: Dirichlet
    - neu: Neumann
    - trans: Transmission

List of wrappers:

- `Helmholtz <helm_wrappers.html#hwrap>`__ 

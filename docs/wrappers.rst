Boundary value problem solvers
===============================

For each of the solver in this repository, we have the following 
user callable routines::

    - lpcomp_<pde>_<rep>_<data>: evaluate layer potential
    - <pde>_<rep>_<bvp>_solver: iterative solver

- <pde>: PDE being solved.

    - helm: Helmholtz
    - lap: Laplace
    - stok: Stokes
    - em: Maxwell

- <rep>: Layer potential representation for solution. Some examples include

    - s: single layer representation
    - comb: combined field representation
    - rpcomb: right preconditioned combied field representation
    - mfie: magnetic field integral representation
    - nrccie: non-resonant charge current integral equation
    - dpie: decoupled potential integral representation
    - muller: adjoint of muller integral equation
    - dfie: decoupled field integral equation

- <bvp/data>: Boundary value problem or post processing data. Some examples include

    - dir: Dirichlet
    - neu: Neumann
    - trans: Acoustic Transmission/Dielectric interface
    - pec: Perfect electric conductor

List of wrappers:

- `Helmholtz <helm_wrappers.html#hwrap>`__ 

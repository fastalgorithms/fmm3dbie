fmm3dBIE MATLAB User Guide
===========================

Apologies that this guide is still under construction. There
are examples for the most common uses. 

Detailed documentation of all high-level functions and
classes can be obtained by viewing the source file or by
typing "help [function name]" in a MATLAB command window:

.. code:: matlab

   help surfer
	  
Currently the following user commands are available::

Defining surface objects from analytical descriptions (+geometries)
--------------------------------------------------------------------

- geometries.sphere
- geometries.ellipsoid
- geomtries.startorus
- geometries.stellarator

Obtaining smooth surface triangulations from low order meshes
--------------------------------------------------------------

- multiscale_mesher

Manipulating surfer objects (@surfer)
------------------------------------------------------

- translate
- rotate
- merge
- affine_transf
- scale

Plotting routines (@surfer)
----------------------------

- plot
- scatter
- plot_nodes

Other useful routines (@surfer)
--------------------------------

- surf_fun_error
- vals2coefs

Solver routines
----------------

- solver

  * lap3d.solver

    - lap3d.dirichlet.solver

    - lap3d.neumann.solver

  * helm3d.solver

    - helm3d.dirichlet.solver

    - helm3d.neumann.solver

    - helm3d.impedance.solver

    - helm3d.transmission.solver

  * stok3d.solver

    - stok3d.velocity.solver

  * em3d.solver

    - em3d.pec.solver


Evaluation routines
--------------------

- eval_fields

  * lap3d.eval

    - lap3d.dirichlet.eval

    - lap3d.neumann.eval

  * helm3d.eval

    - helm3d.dirichlet.eval

    - helm3d.neumann.eval

    - helm3d.impedance.eval

    - helm3d.transmission.eval

  * stok3d.eval

    - stok3d.velocity.eval

  * em3d.eval

    - em3d.pec.eval


Kernel routines
----------------

- lap3d.kern
- helm3d.kern
- helm3d.planewave
- stok3d.kern
- em3d.kern
- em3d.planewave

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   guide/surfers
   guide/solvers
   
   

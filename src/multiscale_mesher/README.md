This repository contains the code for generating high order meshes
for multiscale geometries where the input geometry is defined through a
watertight second order triangulation, and the output is stored
in the go3 format (See
[here](https://fmm3dbie.readthedocs.io/en/latest/surface_reps.html#go3) 
for a description of how the geometry is
stored in the go3 format)

The code currently has no dependencies and the easiest way to run the
code is through the MATLAB app (matlab/Smoother2.mlapp).

If you don't have matlab, you can alternatively run the makefile in the
parent directory. Running the make file runs the test script
`examples/test_surfsmooth.f90`. You can smooth your own mesh, by setting
the string 'nombre' on line 92, and the output file name by setting the
string 'filename'. As in the matlab app, you can also tweak various
parameters of the algorithm above line 92. The output of running the
program is 'nrefine' go3 files whose names will be
`filename_o(#norder_smooth)_r(#refinement number)`. Refinement 0 has the
same number of triangles as the input msh, and every subsequent
refinement has 4 times as many triangles as the previous refinement.

Note that when switching from running the code in
matlab to fortran or vice vrsa, it is essential to run make clean in the
parent directory before the switch.

Tutorials:
* [Smoothing and visualizing multi-component geometries](https://media.upv.es/#/portal/video/a34cf630-7dd8-11eb-8282-fdfd21383b49)

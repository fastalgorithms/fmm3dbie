# fmm3dbie

## FMM-accelerated boundary integral equation solvers

<p align="center">
<img width="50%" src="docs/plane.png"/>
</p>

Currently supports high-order triangulation of smooth surfaces

Upcoming support for: 
-  High order quadrilaterization versions of the above routines 


This repository has an external dependency - [FMM3D](https://fmm3d.readthedocs.io/en/latest)

Make sure you have the shared object for the FMM library installed and
located in an appropriate location (`/usr/local/lib` on MacOSX, and
environment variable of LD_LIBRARY_PATH set to location of libfmm3d.so 
on linux machines)


Please see the [online documentation](https://fmm3dbie.readthedocs.io).


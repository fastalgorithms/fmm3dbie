# fmm3dbie

## FMM-accelerated boundary integral equation solvers

<p align="center">
<img width="50%" src="docs/plane.png"/>
</p>

Currently supports high-order triangulations and quadrilaterizations 
of smooth surfaces

This repository has an external dependency - [FMM3D](https://fmm3d.readthedocs.io/en/latest)

The package along with the dependency can be obtained by


    git clone --recurse-submodules https://github.com/fastalgorithms/fmm3dbie.git


Make sure you have the shared object for the FMM library installed and
located in an appropriate location (`/usr/local/lib` on MacOSX, and
environment variable of LD_LIBRARY_PATH set to location of libfmm3d.so 
on linux machines)


Please see the [online documentation](https://fmm3dbie.readthedocs.io).


fmm3dbie team
===============
* Travis Askham
* Leslie Greengard
* Jeremy Hoskins
* Libin Lu
* Mike O'Neil
* Manas Rachh
* Felipe Vico
* Vladimir Rokhlin


James Bremer provided generalized Gaussian quadrature rules
(`src/quadratures/ggq-self-quadrouts.f`, `src/quadratures/ggq-selfquad.f`,
and all files in `src/quadratures/ggq-self-quads/`). 

Zydrunas Gimbutas provided the high order quadrature rules
for integrating smooth functions `src/tria_routs/triasymq.f`, 
`src/tria_routs/koorn-uvs-dat.txt`, `src/tria_routs/koorn-wts-dat.txt`,
and, `src/quad_routs/squarearbq.f`; and also provided `src/common/dotcross3d.f`

References
============
If you find fmm3dbie useful in your work, please cite this repository,
our paper and the papers associated with quadrature contributions by
James Bremer and Zydrunas Gimbutas:

- `Greengard, L., O'Neil, M., Rachh, M., & Vico, F. (2021). Fast multipole methods for 
the evaluation of layer potentials with locally-corrected quadratures. 
Journal of Computational Physics: X, 10, 100092.`

- `Bremer, J., & Gimbutas, Z. (2013). On the numerical evaluation of the singular 
integrals of scattering theory. Journal of Computational Physics, 251, 327-343.` 

- `Xiao, H., & Gimbutas, Z. (2010). A numerical algorithm for the construction of 
efficient quadrature rules in two and higher dimensions. Computers & mathematics with 
applications, 59(2), 663-676.` 


Funding
=========

Funded in part by the Office of Naval Research under Awards
#N00014-17-1-2059, #N00014-17-1-2451, and #N00014-18-1-2307, 
and the Simons Foundation/SFARI (560651, AB).


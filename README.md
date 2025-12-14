# fmm3dbie

## FMM-accelerated boundary integral equation solvers

<p align="center">
<img width="50%" src="docs/plane.png"/>
</p>

Currently supports high-order triangulations and quadrilaterizations 
of smooth surfaces

This repository has an external dependency - [FMM3D](https://fmm3d.readthedocs.io/en/latest)
which is included as a submodule

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
* Zydrunas Gimbutas
* Tristan Goodwill
* Leslie Greengard
* Jeremy Hoskins
* Libin Lu
* Mike O'Neil
* Manas Rachh
* Vladimir Rokhlin
* Felipe Vico


References
============
If you find fmm3dbie useful in your work, we ask that you please cite the
following works

- {This software} see CITATIONS.cff for details
- {Locally corrected quadratures framework} "Greengard, L., O'Neil, M., Rachh, M., & Vico, F. (2021). Fast multipole methods for 
the evaluation of layer potentials with locally-corrected quadratures. 
Journal of Computational Physics: X, 10, 100092."
- {Quadrature generation routines} "Bremer, J., & Gimbutas, Z. (2013). On the numerical evaluation of the singular 
integrals of scattering theory. Journal of Computational Physics, 251, 327-343." 
- {Quadrature and interpolation nodes for smooth functions on triangles} "Vioreanu, B., & Rokhlin, V. (2014). Spectra of multiplication operators as a numerical tool. SIAM Journal on Scientific Computing 36.1, A267-A288."
- {Quadratures for smooth functions on triangles} "Xiao, H., & Gimbutas, Z. (2010). A numerical algorithm for the construction of 
efficient quadrature rules in two and higher dimensions. Computers & mathematics with 
applications, 59(2), 663-676." 


Funding
=========

Funded in part by the Office of Naval Research under Awards
#N00014-17-1-2059, #N00014-17-1-2451, and #N00014-18-1-2307, 
and the Simons Foundation/SFARI (560651, AB).

Contributing
=============

Contributions are welcome. See the issues tab or create a new issue
if you're interested in bringing something in to fmm3dbie. 
See the [wiki](https://github.com/fastalgorithms/fmm3dbie/wiki)
for more on the developer process.

Acknowledgements
=================
We would also like to acknowledge contributions from 
* Alex Barnett (documentation and examples for the matlab interface) 
* Ludvig af Klinteberg, Peter Nekrasov, and David Krantz (for contributions and bugfixes to the matlab interface)
* Dan Fortunato (for mex compilation)

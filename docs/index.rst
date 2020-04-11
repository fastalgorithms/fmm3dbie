.. FMM3D documentation master file, created by
   sphinx-quickstart on Wed Nov  1 16:19:13 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Boundary integral solvers in 3D (BIE3D)
==================================================

.. image:: plane.png
    :width: 60%
    :align: center
	    
`BIE3D <https://gitlab.com/fastalgorithms/solvers3d>`_ 
is a set of libraries to solve elliptic boundary value
problems on surfaces inthree dimensions.
The library currently supports Dirichlet and Neumann boundary value
problems for Laplace, Helmholtz, Yukawa, Maxwell and Stokes
equations on multi-core shared memory machines.
The library provides support for evaluating layer potentials on and
off the surface using locally corrected precomputed quadratures, 
FMM accelerated iterative solvers, interfaces for matrix entry
generation for dense and fast direct solvers.
The library is written in Fortran,
and has wrappers for C, and Python.
As an example, given a domain $\Omega \in \mathbb{R}^{3}$ 
whose boundary is $\Gamma$,
and a function $g$ defined on $\Gamma$, the interior Dirichlet
problem for Helmholtz's equation is given by

.. math:: 
   
   (\Delta + k^2) u &= 0 \quad \mbox{ in } \Omega \, , \\
   u &= g \quad \mbox{ on } \Gamma \, .

We represent the solution $u$ in $\Omega$ using the combined field
representation 

.. math:: u = ik\mathcal{S}_{k}[\sigma] + \mathcal{D}_{k}[\sigma] \, ,

where $\sigma$ is the unknown density to be solved for, and 
$\mathcal{S}_{k}[\sigma]$ and $\mathcal{D}_{k}[\sigma]$ are the
Helmholtz single and double layer potentials given by

.. math::
   
   \mathcal{S}_{k}[\sigma](x) &= \int_{\Gamma}
   \frac{e^{ik\|x-y\|}}{\|x-y\|} \sigma(y) dS_{y} \\
   \mathcal{D}_{k}[\sigma](x) &= \int_{\Gamma}
   \nabla_{y} \frac{e^{ik\|x-y\|}}{\|x-y\|} \cdot n(y) \sigma(y) dS_{y} \, ,

and $n(y)$ is the normal to the surface $\Gamma$ at $y$. On imposing the
boundary conditions, we get the following integral equation for 
the unknown density $\sigma$

.. math::
   
   -\frac{1}{2} \sigma + ik S_{k}[\sigma] + D_{k}[\sigma] = g \quad
   x \in \Gamma \, .


.. toctree::
   :maxdepth: 2
	   
   install
   

   

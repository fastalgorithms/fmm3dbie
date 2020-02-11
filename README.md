# solvers3d

## FMM-accelerated boundary integral equation solvers


- Currently only supports high-order triangulation of smooth surfaces

Upcoming support for: 
    *  High order triangulation of surfaces with edges
    *  High order triangulation of surfaces with edges and corners
    *  High order quadrilaterization versions of the above routines 


This repository has two external dependencies - FMM3D (enter link here)
and utils (enter link here)

To compile a local archive and shared object for the FMM libraries
run ``make -f make_fmm_lib <options>''  options: "OMP=ON", "FAST_KER=ON"


The surface is assumed to given as a collection of high order
patches $`S = \cup_{j} \Gamma_{j}`$ with a given parameterization

$`xyz_{j}: T_{0} \to \Gamma_{j}`$

srcvals(12,*) - at each discretization node, stores
  x,y,z,dx/du,dy/du,dz/du,dx/dv,dy/dv,dz/dv,nx,ny,nz

srccoefs(9,*) - for each patch store the orthogonal polynomial expansion
   coefficients of x,y,z,dx/du,dy/du,dz/du,dx/dv,dy/dv,dz/dv

ixyzs(npatches+1) - ixyzs(i) indicates the starting point in the srcvals/
srccoefs array for information related to patch i. The number of 
points on the patch/coefs of orthogonal polynomial expansions is 
then given by ixyzs(i+1)-ixyzs(i)


norders(npatches) - order of discretization for each patch

iptype(npatches) - indicates the type of patch
   iptype = 1, triangular patch discretized using Vioreanu Rokhlin nodes


To run example solver for the dirichlet problem for Helmholtz equation,
run the following set of commands
<Enter file name>/Command

For a sample call to various quadarature generation routines and building
fast layer potential matvec routines checkout 

<Enter file name>
  

import numpy as np
import numpy.linalg as la
import fmm3dbie as h3
import fmm3dpy as fmm3d

# setup the geomtry for sphere
# number of refinements
nref = 1

# order of discretization
norder = 3


# The following derived quantites must not be changed
npatches = int(12*4**int(nref))
npols = int((norder+1)*(norder+2)/2)
npts = int(npatches*npols)


norders,ixyzs,iptype,srcvals,srccoefs,wts = h3.get_sphere_geom(nref,
  npatches,norder,npts)

# generate a point in the exterior and exterior
xyz_in = np.random.rand(3)
xyz_in = xyz_in/np.linalg.norm(xyz_in)*0.3

xyz_out = np.random.rand(3)
xyz_out = xyz_out/np.linalg.norm(xyz_out)*2.6

print("error in area of sphere = ",str(sum(wts)-4*np.pi))

sigout = np.random.rand(3)

# flag for deciding whether to solve interior problem or 
# exteior one. ifinout = 0 ==> interior problem
# ifinout = 1 ==> exterior problem

ifinout = 1
if(ifinout == 0):
    xyz_src = xyz_out
    xyz_targ = xyz_in
else:
    xyz_src = xyz_in
    xyz_targ = xyz_out

#
# get boundary data due to a point stokeslet in the appropriate region
#

rhs = np.zeros(3*npts)
# dpars, ipars, and zk are not used in this routine but are dummy variables
ipars = np.zeros(1,dtype='int')
dpars = 0
zk = complex(0.0,0.0)
ndz = 1
np.shape(sigout)
for i in range(npts):
    smat = h3.st3d_slp_vec(9,xyz_src,srcvals[:,i],dpars,ndz,zk,ipars)
    smat = np.reshape(smat,(3,3),'F')
    istart = i*3
    iend = (i+1)*3
    rhs[istart:iend] = np.dot(smat,sigout)


# Set parameter for solver
# alpha is the strength of stokeslet
# beta is the strength of the stresslet
alpha = 1.0
beta = 1.0
dpars = np.array([alpha,beta])
eps = 0.51e-6

# max number of iterations
numit = 200

# gmres residual requirement
eps_gmres = 0.51e-9

# Solve the boundary value problem using iterative solver
niter,errs,res,soln = h3.stok_comb_vel_solver(norders,ixyzs,iptype,srccoefs,
  srcvals,eps,dpars,numit,ifinout,rhs,eps_gmres)

# Compute exact solution
smat = h3.st3d_slp_vec(9,xyz_src,xyz_targ,dpars,ndz,zk,ipars)
smat = np.reshape(smat,(3,3),'F')
uex = np.dot(smat,sigout)

# compute velocity at target using layer potential evaluator routine
#
# ipatch_id is an integer array of size ntarg
# uvs_targ is a double array of size(2,ntarg)
#
# ipatch_id indicated if a target is on the surface or not
# if it is, then only the pv part of the velocity is returned
#
# ipatch_id = -1 otherwise
#
# if target is on surface, then uvs_targ denotes the local uv coordinates
# of a target on a patch, else it is irrelevant

ipatch_id = -1
uvs_targ = np.zeros(2)
soln = np.reshape(soln,(3,npts),'F')
u = h3.lpcomp_stok_comb_vel(norders,ixyzs,iptype,srccoefs,srcvals,xyz_targ,
  ipatch_id,uvs_targ,eps,dpars,soln)

err = np.linalg.norm(np.transpose(u) - uex)
print("error in solution = ",str(err))

# plot magnitude of density on surface
# individual components are also stored in case you want to plot the
# vector field
h3.surf_vtk_plot_vec(norders,ixyzs,iptype,srccoefs,srcvals,soln,'sph-stok-sig.vtk','title')


#
#
#  Now compute the solution using direct inversion
#

xmat = h3.stok_comb_vel_matgen(norders,ixyzs,iptype,srccoefs,srcvals,eps,
  dpars,ifinout)

rhsuse = np.reshape(rhs,3*npts)
soln2 = np.linalg.solve(xmat,rhsuse)

ipatch_id = -1
uvs_targ = np.zeros(2)
soln2 = np.reshape(soln2,(3,npts),'F')
eps = 0.51e-9
u = h3.lpcomp_stok_comb_vel(norders,ixyzs,iptype,srccoefs,srcvals,xyz_targ,
  ipatch_id,uvs_targ,eps,dpars,soln2)

err = np.linalg.norm(np.transpose(u) - uex)
print("error in solution from inverting matrix= ",str(err))



#
# given density evaluate density on the sphere
#
ipatch_id_src,uvs_src = h3.get_patch_id_uvs(norders,ixyzs,iptype,npts)

usurf = h3.lpcomp_stok_comb_vel(norders,ixyzs,iptype,srccoefs,srcvals,
  srcvals,ipatch_id_src,uvs_src,eps,dpars,soln2)

#
# Note routine only returns principal value part of the double
# layer potential on surface. Fix the identity term
#
usurf = usurf - soln2*dpars[1]*(-1.0)**(ifinout)*2*np.pi
rhs = np.reshape(rhs,(3,npts),'F')
err = np.linalg.norm((usurf[0,:]-rhs[0,:])*np.sqrt(wts))/np.linalg.norm(rhs[0,:]*np.sqrt(wts))
print("error in computing velocity on boundary = ",str(err))




import numpy as np
import numpy.linalg as la
import fmm3dbie as h3
import fmm3dpy as fmm3d

x = np.loadtxt('../geometries/sphere_192_o03.go3')

ifwrite = 1

norder = int(x[0])
npatches = int(x[1])
npols = int((norder+1)*(norder+2)/2)
npts = npatches*npols 

# setup geometry in the correct format
norders = norder*np.ones(npatches)
iptype = np.ones(npatches)
srcvals = x[2::].reshape(12,npts)
ixyzs = np.arange(npatches+1)*npols+1

# convert values to coefs
srccoefs = h3.surf_vals_to_coefs(norders,ixyzs,iptype,srcvals[0:9,:])

wts = h3.get_qwts(norders,ixyzs,iptype,srcvals)

print("error in area of sphere = ",str(sum(wts)-4*np.pi))

# using two sources in exterior
# note: code will currently break if one source is used (fmm interace issue)
#
xyz_out = np.array([[3.17,-0.03,3.15],[6.13,-4.1,2.22]]).transpose()
zk = 1.1 + 1j*0
c = np.array([1 + 1j*0,1+1.1j])
out = fmm3d.h3ddir(zk=zk,sources=xyz_out,targets=srcvals[0:3,:],charges=c,pgt=1)
rhs = out.pottarg

alpha = 0.0
beta = -2.0
zpars = np.array([zk,alpha,beta],dtype=complex)
eps = 0.51e-6

nifds,nrfds,nzfds = h3.helm_comb_dir_fds_block_mem(norders,ixyzs,iptype,srccoefs,srcvals,
  eps,zpars,ifwrite)

ifds,zfds = h3.helm_comb_dir_fds_block_init(norders,ixyzs,iptype,srccoefs,srcvals,eps,
  zpars,nifds,nzfds)

row_ind = np.arange(npts)+1
col_ind = np.arange(npts)+1

xmat = h3.helm_comb_dir_fds_block_matgen(norders,ixyzs,iptype,srccoefs,srcvals,
  wts,eps,zpars,ifds,zfds,row_ind,col_ind,ifwrite)

xmat = xmat - np.identity(npts)*zpars[2]/2.0
sigma = la.solve(xmat,rhs)

xyz_in = np.array([[0.17,-0.03,0.15],[0.13,-0.1,0.22]]).transpose()
ntarg = np.shape(xyz_in)[1]

ipatch_id = -1*np.ones(2)
uvs_targ = np.zeros((2,ntarg)) 

pot_comp = h3.lpcomp_helm_comb_dir(norders,ixyzs,iptype,srccoefs,srcvals,
   xyz_in,ipatch_id,uvs_targ,eps,zpars,sigma)

out = fmm3d.h3ddir(zk=zk,sources=xyz_out,targets=xyz_in,charges=c,pgt=1)
pot_ex = out.pottarg
erra = la.norm(pot_ex-pot_comp)
print("error in solution = ",str(erra))

# plot surface with shading as z coordinate
h3.surf_vtk_plot(norders,ixyzs,iptype,srccoefs,srcvals,'sph.vtk','title1')

rsigma = np.real(sigma)
# plot surface with real part of density as solution
h3.surf_vtk_plot_scalar(norders,ixyzs,iptype,srccoefs,srcvals,rsigma,'sph-sig.vtk','title1')


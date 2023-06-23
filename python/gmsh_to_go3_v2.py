# -----------------------------------------------------------------------------
#
#  Gmsh to Go3
#
# -----------------------------------------------------------------------------

# Usage:
# python3 gmsh_to_go3_v2.py gmsh_filename
# output: gmsh_filename.go3

# one needs the kexp and srout module
# run `make python-gmsh` in fmm3dbie root directory

# test
# run `python3 gmsh_to_go3_v2.py gmsh_test1.msh`
# check output gmsh_test1.msh.go3

# The Python script convert gmsh to go3
# supports 
# 1. first order triangle elemtype 2, 3 nodes
# 2. second order triangle elemtype 9, 6 nodes
# 3. third order triangle elemtype 21, 10 nodes
# 4. forth order triangle elemtype 23, 15 nodes
# 5. fifth order triangle elemtype 25, 21 nodes

import gmsh
import sys
import numpy as np
import kexp
import srout
import os

if len(sys.argv) < 2:
    print("Usage: " + sys.argv[0] + " file")
    exit

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

fname = sys.argv[1]+".go3"
try:
    os.remove(fname)
except OSError:
    pass
gmsh.open(sys.argv[1])

# Print the model name and dimension:
print('Model ' + gmsh.model.getCurrent() + ' (' +
      str(gmsh.model.getDimension()) + 'D)')

supportedTypes = [2,9,21,23,25]
ndim = 3

elem_types = gmsh.model.mesh.getElementTypes()
# * List all types of elements making up the mesh of the entity:
for t in elem_types:
    if t not in supportedTypes:
        continue

    name, dim, order, numv, uvs, _ = gmsh.model.mesh.getElementProperties(t)

    print(" - Element type: " + name + ", order " + str(order) + " (" +
          str(numv) + " nodes in param coord: " + str(uvs) + ")")

    elemTags, elemNodeTags = gmsh.model.mesh.getElementsByType(t)

    nume = len(elemTags)
    npols = int((order+1)*(order+2)/2)
    uvs = uvs.reshape([npols,2]).transpose()
    eCoef = np.zeros([ndim, npols*nume])
    ePos = np.zeros([ndim,npols*nume])
    eDu = np.zeros([ndim,npols*nume])
    eDv = np.zeros([ndim,npols*nume])
    eNormal = np.zeros([ndim,npols*nume])
    umatr = kexp.koorn_vals2coefs(nmax=order, npols=npols, uvs=uvs)
    vmatr = kexp.koorn_coefs2vals_vioreanu(norder=order, npols=npols)

    print("processing element type: "+name)
    for ielem in range(0,nume):
        for ipol in range(0,npols):
            for idim in range(0,ndim):
                indi = ipol + ielem*npols
                for jpol in range(0,npols):
                    indj = int(elemNodeTags[jpol+ielem*numv])
                    eCoef[idim][indi] = eCoef[idim][indi] + umatr[ipol][jpol]*gmsh.model.mesh.getNode(indj)[0][idim]

    print("coefs to vals")
    for ielem in range(0,nume):
        for ipol in range(0,npols):
            for idim in range(0,ndim):
                indi = ipol + ielem*npols
                for jpol in range(0,npols):
                    ePos[idim][indi] = ePos[idim][indi] + vmatr[ipol][jpol]*eCoef[idim][jpol+ielem*npols]

    print("compute du, dv and normals")
    ptcnt=0
    for ielem in range(0,nume):
        ex = ePos[0][ielem*npols:(ielem+1)*npols]
        ey = ePos[1][ielem*npols:(ielem+1)*npols]
        ez = ePos[2][ielem*npols:(ielem+1)*npols]
        ex = ex.reshape([1,npols])
        ey = ey.reshape([1,npols])
        ez = ez.reshape([1,npols])
        exduv = srout.get_surf_uv_grad_tri(nd=1,norder=order,npols=npols,f=ex).reshape([2,npols])
        eyduv = srout.get_surf_uv_grad_tri(nd=1,norder=order,npols=npols,f=ey).reshape([2,npols])
        ezduv = srout.get_surf_uv_grad_tri(nd=1,norder=order,npols=npols,f=ez).reshape([2,npols])
        for ipts in range(0,npols):
            dxyzu = [exduv[0][ipts], eyduv[0][ipts], ezduv[0][ipts]]
            dxyzv = [exduv[1][ipts], eyduv[1][ipts], ezduv[1][ipts]]
            normal = np.cross(dxyzu,dxyzv)
            normal = normal/np.sqrt(np.sum(normal*normal))
            eDu[0][ptcnt] = exduv[0][ipts]
            eDu[1][ptcnt] = eyduv[0][ipts]
            eDu[2][ptcnt] = ezduv[0][ipts]
            eDv[0][ptcnt] = exduv[1][ipts]
            eDv[1][ptcnt] = eyduv[1][ipts]
            eDv[2][ptcnt] = ezduv[1][ipts]
            eNormal[0][ptcnt] = normal[0];
            eNormal[1][ptcnt] = normal[1];
            eNormal[2][ptcnt] = normal[2];
            ptcnt=ptcnt+1

    f=open(fname,'a')
    f.write(str(order))
    f.write('\n')
    f.close()
    f=open(fname,'a')
    f.write(str(nume))
    f.write('\n')
    np.savetxt(f, ePos.ravel(), delimiter='\t')
    np.savetxt(f, eDu.ravel(), delimiter='\t')
    np.savetxt(f, eDv.ravel(), delimiter='\t')
    np.savetxt(f, eNormal.ravel(), delimiter='\t')
    f.close()

# We can use this to clear all the model data:
gmsh.clear()

gmsh.finalize()

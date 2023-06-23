# -----------------------------------------------------------------------------
#
#  Gmsh to Go3
#
# -----------------------------------------------------------------------------

# Usage:
# python3 gmsh_to_go3_numba.py gmsh_filename
# output: gmsh_filename.go3

# one needs the kexp and srout module
# run `make python-gmsh` in fmm3dbie root directory

# test
# run `python3 gmsh_to_go3_numba.py gmsh_test1.msh`
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
import numba as nb
import kexp
import srout
import os
import ctypes
from ctypes import pythonapi
from numba import njit, prange,jit
from numba import types
from numba.extending import intrinsic
from numba.core import cgutils

@intrinsic
def val_to_ptr(typingctx, data):
    def impl(context, builder, signature, args):
        ptr = cgutils.alloca_once_value(builder,args[0])
        return ptr
    sig = types.CPointer(nb.typeof(data).instance_type)(nb.typeof(data).instance_type)
    return sig, impl

def capsule_name(capsule):
    pythonapi.PyCapsule_GetName.restype = ctypes.c_char_p
    pythonapi.PyCapsule_GetName.argtypes = [ctypes.py_object]
    return pythonapi.PyCapsule_GetName(capsule)

def get_f2py_function_address(capsule):
    name = capsule_name(capsule)
    pythonapi.PyCapsule_GetPointer.restype = ctypes.c_void_p
    pythonapi.PyCapsule_GetPointer.argtypes = [ctypes.py_object, ctypes.c_char_p]
    return pythonapi.PyCapsule_GetPointer(capsule, name)

dbl_p=ctypes.POINTER(ctypes.c_double)
int_p =ctypes.POINTER(ctypes.c_int)
functype = ctypes.CFUNCTYPE(ctypes.c_void_p,
                            int_p,int_p,int_p,dbl_p,dbl_p)
func_ptr=get_f2py_function_address(srout.get_surf_uv_grad_tri._cpointer)
srout_normal_du_dv = functype(func_ptr)

@njit(parallel=True)
def val2coef(nume,npols,ndim,numv,elemNodeTags,eCoef,umatr,coord):
    for ielem in prange(nume):
        for ipol in range(0,npols):
            for idim in range(0,ndim):
                indi = ipol + ielem*npols
                for jpol in range(0,npols):
                    indj = int(elemNodeTags[jpol+ielem*numv])
                    eCoef[idim][indi] = eCoef[idim][indi] + umatr[ipol][jpol]*coord[idim+(indj-1)*ndim]

@njit(parallel=True)
def coef2val(nume,npols,ndim,ePos,eCoef,vmatr):
    for ielem in prange(nume):
        for ipol in range(0,npols):
            for idim in range(0,ndim):
                indi = ipol + ielem*npols
                for jpol in range(0,npols):
                    ePos[idim][indi] = ePos[idim][indi] + vmatr[ipol][jpol]*eCoef[idim][jpol+ielem*npols]

@njit(parallel=True)
def get_normal_du_dv(nume,npols,order,ePos,eDu,eDv,eNormal):
    nd_ptr = val_to_ptr(nb.int32(1))
    order_ptr = val_to_ptr(nb.int32(order))
    npols_ptr = val_to_ptr(nb.int32(npols))
    for ielem in prange(0,nume):
        ptcnt = ielem*npols
        ex = ePos[0][ielem*npols:(ielem+1)*npols]
        ey = ePos[1][ielem*npols:(ielem+1)*npols]
        ez = ePos[2][ielem*npols:(ielem+1)*npols]
        exduv = np.empty((npols,2),dtype=np.float64)
        eyduv = np.empty((npols,2),dtype=np.float64)
        ezduv = np.empty((npols,2),dtype=np.float64)
        srout_normal_du_dv(nd_ptr,order_ptr,npols_ptr,ex.ctypes,exduv.ctypes)
        srout_normal_du_dv(nd_ptr,order_ptr,npols_ptr,ey.ctypes,eyduv.ctypes)
        srout_normal_du_dv(nd_ptr,order_ptr,npols_ptr,ez.ctypes,ezduv.ctypes)
        for ipts in range(0,npols):
            dxyzu = [exduv[ipts][0], eyduv[ipts][0], ezduv[ipts][0]]
            dxyzv = [exduv[ipts][1], eyduv[ipts][1], ezduv[ipts][1]]
            normal = np.cross(dxyzu,dxyzv)
            normal = normal/np.sqrt(np.sum(normal*normal))
            eDu[0][ptcnt] = exduv[ipts][0]
            eDu[1][ptcnt] = eyduv[ipts][0]
            eDu[2][ptcnt] = ezduv[ipts][0]
            eDv[0][ptcnt] = exduv[ipts][1]
            eDv[1][ptcnt] = eyduv[ipts][1]
            eDv[2][ptcnt] = ezduv[ipts][1]
            eNormal[0][ptcnt] = normal[0];
            eNormal[1][ptcnt] = normal[1];
            eNormal[2][ptcnt] = normal[2];
            ptcnt=ptcnt+1



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
nodes, coord, _ = gmsh.model.mesh.getNodes()

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
    umatr = kexp.koorn_vals2coefs(nmax=order, npols=npols, uvs=uvs)
    vmatr = kexp.koorn_coefs2vals_vioreanu(norder=order, npols=npols)
    eDu = np.zeros([ndim,npols*nume])
    eDv = np.zeros([ndim,npols*nume])
    eNormal = np.zeros([ndim,npols*nume])

    print("processing element type: "+name)
    print("gmsh val to koorn coef: "+name)
    val2coef(nume,npols,ndim,numv,elemNodeTags,eCoef,umatr,coord)
    print("done gmsh val to koorn coef: "+name)

    print("koorn coefs to vals on koorn points")
    coef2val(nume,npols,ndim,ePos,eCoef,vmatr)
    print("done koorn coefs to vals on koorn points")

    print("compute du, dv and normals")
    get_normal_du_dv(nume,npols,order,ePos,eDu,eDv,eNormal)
    print("done compute du, dv and normals")

    print("writing .go3 file")
    f=open(fname,'a')
    f.write(str(order))
    f.write('\n')
    f.close()
    f=open(fname,'a')
    f.write(str(nume))
    f.write('\n')
    np.savetxt(f, ePos.ravel())
    np.savetxt(f, eDu.ravel())
    np.savetxt(f, eDv.ravel())
    np.savetxt(f, eNormal.ravel())
    f.close()
    print("done writing .go3 file")

# We can use this to clear all the model data:
gmsh.clear()

gmsh.finalize()

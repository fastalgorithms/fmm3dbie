# -----------------------------------------------------------------------------
#
#  Gmsh to Go3
#
# -----------------------------------------------------------------------------

# Usage:
# python3 gmsh_to_go3.py gmsh_filename
# output: gmsh_filename.go3

# one needs the kexp and srout module
# run `make python-gmsh` in fmm3dbie root directory

# test
# run `python3 gmsh_to_go3.pyt gmsh_test1.msh`
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

# Geometrical data is made of elementary model `entities', called `points'
# (entities of dimension 0), `curves' (entities of dimension 1), `surfaces'
# (entities of dimension 2) and `volumes' (entities of dimension 3). As we have
# seen in the other Python tutorials, elementary model entities are identified
# by their dimension and by a `tag': a strictly positive identification
# number. Model entities can be either CAD entities (from the built-in `geo'
# kernel or from the OpenCASCADE `occ' kernel) or `discrete' entities (defined
# by a mesh). `Physical groups' are collections of model entities and are also
# identified by their dimension and by a tag.

# Get all the elementary entities in the model, as a vector of (dimension, tag)
# pairs:
entities = gmsh.model.getEntities()

supportedTypes = [2,9,21,23,25]
ndim = 3

for e in entities:
    # Dimension and tag of the entity:
    dim = e[0]
    tag = e[1]

    # Mesh data is made of `elements' (points, lines, triangles, ...), defined
    # by an ordered list of their `nodes'. Elements and nodes are identified by
    # `tags' as well (strictly positive identification numbers), and are stored
    # ("classified") in the model entity they discretize. Tags for elements and
    # nodes are globally unique (and not only per dimension, like entities).

    # A model entity of dimension 0 (a geometrical point) will contain a mesh
    # element of type point, as well as a mesh node. A model curve will contain
    # line elements as well as its interior nodes, while its boundary nodes will
    # be stored in the bounding model points. A model surface will contain
    # triangular and/or quadrangular elements and all the nodes not classified
    # on its boundary or on its embedded entities. A model volume will contain
    # tetrahedra, hexahedra, etc. and all the nodes not classified on its
    # boundary or on its embedded entities.

    # Get the mesh nodes for the entity (dim, tag):
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim, tag)

    # Get the mesh elements for the entity (dim, tag):
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)

    # Elements can also be obtained by type, by using `getElementTypes()'
    # followed by `getElementsByType()'.

    # Let's print a summary of the information available on the entity and its
    # mesh.

    # * Type and name of the entity:
    typeE = gmsh.model.getType(e[0], e[1])
    name = gmsh.model.getEntityName(e[0], e[1])
    if len(name): name += ' '
    print("Entity " + name + str(e) + " of type " + typeE)

    # * Number of mesh nodes and elements:
    numElem = sum(len(i) for i in elemTags)
    print(" - Mesh has " + str(len(nodeTags)) + " nodes and " + str(numElem) +
          " elements")

    # * Entities on its boundary:
    boundary = gmsh.model.getBoundary([e])
    if len(boundary):
        print(" - Boundary entities: " + str(boundary))

    # * Does the entity belong to physical groups?
    physicalTags = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
    if len(physicalTags):
        s = ''
        for p in physicalTags:
            n = gmsh.model.getPhysicalName(dim, p)
            if n: n += ' '
            s += n + '(' + str(dim) + ', ' + str(p) + ') '
        print(" - Physical groups: " + s)

    # * Is the entity a partition entity? If so, what is its parent entity?
    partitions = gmsh.model.getPartitions(e[0], e[1])
    if len(partitions):
        print(" - Partition tags: " + str(partitions) + " - parent entity " +
              str(gmsh.model.getParent(e[0], e[1])))

    # * List all types of elements making up the mesh of the entity:
    elemTypeIter = 0
    for t in elemTypes:
        if t not in supportedTypes:
            continue
        print(t == elemTypes[elemTypeIter])
        name, dim, order, numv, uvs, _ = gmsh.model.mesh.getElementProperties(
            t)
        print(" - Element type: " + name + ", order " + str(order) + " (" +
              str(numv) + " nodes in param coord: " + str(uvs) + ")")

        nume = len(elemTags[elemTypeIter])
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
                        indj = int(elemNodeTags[elemTypeIter][jpol+ielem*numv])
                        eCoef[idim][indi] = eCoef[idim][indi] + umatr[ipol][jpol]*nodeCoords[idim+(indj-1)*ndim]

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



#        print("debug compute du, dv coefs")
#        print(elemNodeTags[elemTypeIter])
#        print(elemTags[elemTypeIter])
#        print(elemTypes[elemTypeIter])
#        print("end of processing element type: "+name)

        elemTypeIter = elemTypeIter + 1
        f=open(fname,'x')
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

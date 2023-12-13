
//
// This file contains the routines for computing a "skeleton mesh"
// from a "base mesh". The skeleton mesh is basically just the base
// mesh except it has been discretized with Rokhlin-Vioreanu nodes,
// and possibly over-sampled.
//
//



#include "cisurf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cprini.h"



void makeSkeleton( baseMesh *baseMesh1, skelMesh *skelMesh1, long norder ) {
  // constructs a skeleton mesh from the base mesh by popping down
  // Vioreanu-Rokhlin quadrature nodes

  double *uvs;
  double *whts;
  long npols;
  
  // get the reference quadrature notes and weights
  
  // TODO this should probably be either in C, or have some proper
  // c-datatype bindings
  npols = (norder+1)*(norder+2)/2;
  uvs = (double *) malloc( 2*npols*sizeof(double) );
  whts = (double *) malloc( npols*sizeof(double) );
  get_vioreanu_nodes_wts_( &norder, &npols, uvs, whts );

  cprind_matrix("Vioreanu-Rokhlin nodes =", uvs, npols, 2);

  // now cycle through all each baseMesh element and put nodes on it
  long nelems = baseMesh1->nelems;
  cprinf("number of elements = ", &nelems, 1);

  long i, j, nv, iv;
  for (i=0; i<nelems; i++) {
    // print out the info for this element to check things
    printBaseElementInfo( &(baseMesh1->elements[i]) );

    nv = baseMesh1->elements[i].nv;
    cprinf("nv = ", &nv, 1);

    for (j=0; j<nv; j++) {
      iv = baseMesh1->elements[i].ivs[j];
      cprind("vert = ", &(baseMesh1->verts[3*iv]), 3);
    }

    exit(0);
  }


  free( uvs );
  free( whts );

  return;'
}





void eval_base_element( double *uv, baseElement *element, pointInfo *xyz_info) {
  // evaluate the base element at the point uv, proceed by figuring out the type
  // of base element and then going down and calling the guru routine and
  // finally packing everything up


}





void eval_tria2( double *uv, double *verts, double *xyz, double *du, double *dv,
                 double *normal ) {
  // eval the tria2 element using the vertex values

  // use direct formula


  return;
}

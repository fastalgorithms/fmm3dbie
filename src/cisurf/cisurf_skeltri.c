
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
  int npols;
  
  // get the reference quadrature notes and weights
  
  // TODO this should probably be either in C, or have some proper
  // c-datatype bindings
  npols = (norder+1)*(norder+2)/2;
  uvs = (double *) malloc( 2*npols*sizeof(double) );
  whts = (double *) malloc( npols*sizeof(double) );
  get_vioreanu_nodes_wts_( &norder, &npols, uvs, whts );

  cprind_matrix("Vioreanu-Rokhlin nodes =", uvs, npols, 2);

  // now cycle through all each baseMesh element and put nodes on it




  free( uvs );
  free( whts );

}






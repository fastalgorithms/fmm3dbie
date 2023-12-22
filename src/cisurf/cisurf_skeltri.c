
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
#include "cprini_long.h"



void create_skeleton( BaseMesh *basemesh1, SkelMesh *skelmesh1, long norder ) {
  // constructs a skeleton mesh from the base mesh by popping down
  // Vioreanu-Rokhlin quadrature nodes

  double *uvs;
  double *whts;
  long npols;
  
  // get the reference quadrature notes and weights
  //
  // TODO this should probably be either in C, or have some proper
  // c-datatype bindings

  npols = (norder+1)*(norder+2)/2;
  uvs = (double *) malloc( 2*npols*sizeof(double) );
  whts = (double *) malloc( npols*sizeof(double) );
  get_vioreanu_nodes_wts_( &norder, &npols, uvs, whts );

  cprind_matrix("Vioreanu-Rokhlin nodes =", uvs, npols, 2);

  long nelems = basemesh1->nelems;
  cprinf("number of elements = ", &nelems, 1);


  // we'll create a skeleton mesh with THE SAME number of elements on it for
  // now, refinement is done after the fact
  skelmesh1->nelems = nelems;
  skelmesh1->elements = (SkelElement *) malloc( nelems*sizeof(SkelElement) );

  // create a shortcut pointer for these elements
  SkelElement *elements;
  elements = skelmesh1->elements;


  long i, j, nv, iv, nuv, k, l;
  double uv[2], verts[100], xyz[3], du[3], dv[3], da, normal[3];

  double xyzs[1000], dus[1000], dvs[1000];
  double das[1000], normals[1000];

  // now cycle through all each basemesh element and put nodes on it
  for (i=0; i<nelems; i++) {
    // print out the info for this element to check things
    //print_base_element_info( &(basemesh1->elements[i]) );

    nv = basemesh1->elements[i].nv;
    //cprinf("nv = ", &nv, 1);

    for (j=0; j<nv; j++) {
      iv = basemesh1->elements[i].ivs[j];
      verts[3*j] = basemesh1->verts[3*iv];
      verts[3*j+1] = basemesh1->verts[3*iv+1];
      verts[3*j+2] = basemesh1->verts[3*iv+2];
    }

    //cprind_matrix("all the vertices = ", verts, nv, 3);

    // TODO put in switch statement for various kinds of base elements, tria,
    // quad, etc.

    eval_tria2( npols, uvs, verts, xyzs, dus, dvs, das, normals);

    //cprind_matrix("after eval_tria2, xyzs = ", xyzs, npols, 3);
    //cprind_matrix("after eval_tria2, dus = ", xyzs, npols, 3);
    //cprind_matrix("after eval_tria2, dvs = ", xyzs, npols, 3);
    //cprind("after eval_tria2, area elements are = ", das, npols);
    //cprind_matrix("after eval_tria2, normals = ", xyzs, npols, 3);

    //cprin_skipline(2);

    // now copy this stuff over to a skeleton mesh element
    elements[i].id = basemesh1->elements[i].id;
    elements[i].gtype = "tria2";
    elements[i].nv = nv;
    elements[i].verts = (double *) malloc( 3*nv*sizeof(double) );

    for (j=0; j<(3*nv); j++) {
      elements[i].verts[j] = verts[j];
    }

    cprind("all the verts = ", verts, 3*nv);


    // copy of some quadrature weight
    elements[i].norder = norder;
    elements[i].npols = npols;
    elements[i].whts = (double *) malloc(npols*sizeof(double));

    for (j=0; j<npols; j++) {
      elements[i].whts[j] = das[j]*whts[j];
    }

    // copy over some quadrature node and geometry information
    elements[i].srcvals = (PointInfo *) malloc( npols*sizeof(PointInfo) );

    for (j=0; j<npols; j++) {
      for (k=0; k<3; k++) {
        elements[i].srcvals[j].xyz[k] = xyzs[3*j+k];
        elements[i].srcvals[j].du[k] = dus[3*j+k];
        elements[i].srcvals[j].dv[k] = dvs[3*j+k];
        elements[i].srcvals[j].normal[k] = normals[3*j+k];
      }
    }

    // don't eval the coefs yet...
    elements[i].coefs = NULL;

    // some additional geometry info
    elements[i].centroid[0] = basemesh1->elements[i].centroid[0];
    elements[i].centroid[1] = basemesh1->elements[i].centroid[1];
    elements[i].centroid[2] = basemesh1->elements[i].centroid[2];
    elements[i].radius = basemesh1->elements[i].radius;

    // print out this element information to make sure we're loading it properly
    //print_skeleton_element_info( &elements[i] );
    //

  }




  free( uvs );
  free( whts );

  return;
}





void print_skeleton_element_info( SkelElement *element ) {

  // print out the data for a skeleton element

  printf("\n");
  printf("- - - - skel element information - - - - \n");
  printf("element.id          = %ld\n", element->id);
  printf("element.gtype       = %s\n", element->gtype);
  printf("element.nv          = %ld\n", element->nv);

  long i;
  for (i=0; i<element->nv; i++){
    printf("element.verts[%ld]    = (%e, %e, %e)\n", i, element->verts[3*i],
           element->verts[3*i+1], element->verts[3*i+2]);
  }

  printf("element.norder      = %ld\n", element->norder);
  printf("element.npols       = %ld\n", element->npols);

  // don't print out the nodes and weights just yet...

  printf("element.centroid    = (%e, %e, %e)\n", element->centroid[0],
         element->centroid[1], element->centroid[2]);
  printf("element.radius      = %e\n", element->radius);
  printf("- - - end skel element information - - - \n");
  printf("\n");

}





void eval_base_element( double *uv, BaseElement *element, PointInfo *xyzInfo) {
  // evaluate the base element at the point uv, proceed by figuring out the type
  // of base element and then going down and calling the guru routine and
  // finally packing everything up


}





void eval_tria2( long n, double *uv, double *verts, double *xyz, double *du,
                 double *dv, double *da, double *normal ) {
  // eval the tria2 element using the vertex locations
  //
  // Input:
  //   uv - the point on the simplex at which to evaluate the information
  //   verts - input vertices, they are assumed to be ordered as:
  //                 3
  //                 |  \
  //                 6    5
  //                 |       \
  //                 1--- 4 ---2
  //
  // Output:
  //   xyz - the point in R^3 on the quadratic triangle
  //   du - dxdu, dydu, and dzdu
  //   dv - dxdv, dydv, and dzdv
  //   normal - the unit normal vector
  //
  // use direct formula
  //
  long i;
  double coefs[6];
  double p1, p2, p3, p4, p5, p6, u, v;

  for (i=0; i<n; i++) {

    u = uv[2*i];
    v = uv[2*i+1];

    // evaluate the x values first
    p1 = verts[0];
    p2 = verts[3];
    p3 = verts[6];
    p4 = verts[9];
    p5 = verts[12];
    p6 = verts[15];

    coefs[0] = p1;
    coefs[1] = -3*p1 - p2 + 4*p4;
    coefs[2] = -3*p1 - p3 + 4*p6;
    coefs[3] = 2*p1 + 2*p2 - 4*p4;
    coefs[4] = 2*p1 + 2*p3 - 4*p6;
    coefs[5] = 4*p1 - 4*p4 + 4*p5 - 4*p6;

    xyz[3*i] = coefs[0] + coefs[1]*u + coefs[2]*v + coefs[3]*u*u
      + coefs[4]*v*v + coefs[5]*u*v;

    du[3*i] = coefs[1] + 2*coefs[3]*u + coefs[5]*v;
    dv[3*i] = coefs[2] + 2*coefs[4]*v + coefs[5]*u;


    // and now the y values
    p1 = verts[1];
    p2 = verts[4];
    p3 = verts[7];
    p4 = verts[10];
    p5 = verts[13];
    p6 = verts[16];

    coefs[0] = p1;
    coefs[1] = -3*p1 - p2 + 4*p4;
    coefs[2] = -3*p1 - p3 + 4*p6;
    coefs[3] = 2*p1 + 2*p2 - 4*p4;
    coefs[4] = 2*p1 + 2*p3 - 4*p6;
    coefs[5] = 4*p1 - 4*p4 + 4*p5 - 4*p6;

    xyz[3*i+1] = coefs[0] + coefs[1]*u + coefs[2]*v + coefs[3]*u*u
      + coefs[4]*v*v + coefs[5]*u*v;

    du[3*i+1] = coefs[1] + 2*coefs[3]*u + coefs[5]*v;
    dv[3*i+1] = coefs[2] + 2*coefs[4]*v + coefs[5]*u;

    // and now the z values
    p1 = verts[2];
    p2 = verts[5];
    p3 = verts[8];
    p4 = verts[11];
    p5 = verts[14];
    p6 = verts[17];

    coefs[0] = p1;
    coefs[1] = -3*p1 - p2 + 4*p4;
    coefs[2] = -3*p1 - p3 + 4*p6;
    coefs[3] = 2*p1 + 2*p2 - 4*p4;
    coefs[4] = 2*p1 + 2*p3 - 4*p6;
    coefs[5] = 4*p1 - 4*p4 + 4*p5 - 4*p6;

    xyz[3*i+2] = coefs[0] + coefs[1]*u + coefs[2]*v + coefs[3]*u*u
      + coefs[4]*v*v + coefs[5]*u*v;

    du[3*i+2] = coefs[1] + 2*coefs[3]*u + coefs[5]*v;
    dv[3*i+2] = coefs[2] + 2*coefs[4]*v + coefs[5]*u;

    // and finally compute the normal, and the normalize it
    cross_product_3d(&du[3*i], &dv[3*i], &normal[3*i]);

    da[i] = sqrt( normal[3*i]*normal[3*i] + normal[3*i+1]*normal[3*i+1]
                  + normal[3*i+2]*normal[3*i+2] );
    normal[3*i] = normal[3*i]/(da[i]);
    normal[3*i+1] = normal[3*i+1]/(da[i]);
    normal[3*i+2] = normal[3*i+2]/(da[i]);

  }

  return;

}






void eval_tria1( long n, double *uv, double *verts, double *xyz, double *du,
                 double *dv, double *da, double *normal ) {
  // eval the tria1 element using the vertex locations
  //
  // Input:
  //   uv - the point on the simplex at which to evaluate the information
  //   verts - input vertices, they are assumed to be ordered as:
  //                 3
  //                 |  \
  //                 |    \
  //                 |       \
  //                 1---------2
  //
  // Output:
  //   xyz - the point in R^3 on the quadratic triangle
  //   du - dxdu, dydu, and dzdu
  //   dv - dxdv, dydv, and dzdv
  //   normal - the unit normal vector
  //
  // use direct formula
  //

  return;
}





void cross_product_3d(double *x, double *y, double *z) {
  // this routine computes z = x X y, where X is the cross product

  z[0] = x[1]*y[2] - x[2]*y[3];
  z[1] = x[2]*y[0] - x[0]*y[3];
  z[2] = x[0]*y[1] - x[1]*y[0];

  return;
}

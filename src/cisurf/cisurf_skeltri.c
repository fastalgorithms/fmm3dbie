
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

  long i, j, nv, iv, nuv, k, l;
  double uv[2], verts[18], xyz[3], du[3], dv[3], da, normal[3];

  for (i=0; i<nelems; i++) {
    // print out the info for this element to check things
    printBaseElementInfo( &(baseMesh1->elements[i]) );

    nv = baseMesh1->elements[i].nv;
    cprinf("nv = ", &nv, 1);

    for (j=0; j<nv; j++) {
      iv = baseMesh1->elements[i].ivs[j];
      cprind("vert = ", &(baseMesh1->verts[3*iv]), 3);





    }

    cprin_skipline(2);

    // load up the verts
    for (k=0; k<6; k++) {
      iv = baseMesh1->elements[i].ivs[k];
      verts[3*k] = baseMesh1->verts[3*iv];
      verts[3*k+1] = baseMesh1->verts[3*iv+1];
      verts[3*k+2] = baseMesh1->verts[3*iv+2];
    }

    nuv = 1;

    // eval the triangle
    uv[0] = 0;
    uv[1] = 0;
    cprind("uv point = ", uv, 2);
    eval_tria2( nuv, uv, verts, xyz, du, dv, da, normal);
    cprind("from eval_tria2, xyz = ", xyz, 3);

    uv[0] = 1;
    uv[1] = 0;
    cprind("uv point = ", uv, 2);
    eval_tria2( nuv, uv, verts, xyz, du, dv, da, normal);
    cprind("from eval_tria2, xyz = ", xyz, 3);

    uv[0] = 0;
    uv[1] = 1;
    cprind("uv point = ", uv, 2);
    eval_tria2( nuv, uv, verts, xyz, du, dv, da, normal);
    cprind("from eval_tria2, xyz = ", xyz, 3);

    uv[0] = 0.5;
    uv[1] = 0;
    cprind("uv point = ", uv, 2);
    eval_tria2( nuv, uv, verts, xyz, du, dv, da, normal);
    cprind("from eval_tria2, xyz = ", xyz, 3);

    uv[0] = 0.5;
    uv[1] = 0.5;
    cprind("uv point = ", uv, 2);
    eval_tria2( nuv, uv, verts, xyz, du, dv, da, normal);
    cprind("from eval_tria2, xyz = ", xyz, 3);

    uv[0] = 0;
    uv[1] = 0.5;
    cprind("uv point = ", uv, 2);
    eval_tria2( nuv, uv, verts, xyz, du, dv, da, normal);
    cprind("from eval_tria2, xyz = ", xyz, 3);

    exit(0);
  }


  free( uvs );
  free( whts );

  return;
}





void eval_base_element( double *uv, baseElement *element, pointInfo *xyz_info) {
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
  double coefs[6];
  double p1, p2, p3, p4, p5, p6, u, v;

  u = uv[0];
  v = uv[1];

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

  xyz[0] = coefs[0] + coefs[1]*u + coefs[2]*v + coefs[3]*u*u
    + coefs[4]*v*v + coefs[5]*u*v;

  du[0] = coefs[1] + 2*coefs[3]*u + coefs[5]*v;
  dv[0] = coefs[2] + 2*coefs[4]*v + coefs[5]*u;


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

  xyz[1] = coefs[0] + coefs[1]*u + coefs[2]*v + coefs[3]*u*u
    + coefs[4]*v*v + coefs[5]*u*v;

  du[1] = coefs[1] + 2*coefs[3]*u + coefs[5]*v;
  dv[1] = coefs[2] + 2*coefs[4]*v + coefs[5]*u;

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

  xyz[2] = coefs[0] + coefs[1]*u + coefs[2]*v + coefs[3]*u*u
    + coefs[4]*v*v + coefs[5]*u*v;

  du[2] = coefs[1] + 2*coefs[3]*u + coefs[5]*v;
  dv[2] = coefs[2] + 2*coefs[4]*v + coefs[5]*u;

  // and finally compute the normal, and the normalize it
  cross_product_3d(du, dv, normal);

  *da = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
  normal[0] = normal[0]/(*ds);
  normal[1] = normal[1]/(*ds);
  normal[2] = normal[2]/(*ds);

  return;

}

/*
    coef_y(1)=P1(2)
    coef_y(2)=-3*P1(2)-P2(2)+4*P4(2)
    coef_y(3)=-3*P1(2)-P3(2)+4*P6(2)
    coef_y(4)=2*P1(2)+2*P2(2)-4*P4(2)
    coef_y(5)=2*P1(2)+2*P3(2)-4*P6(2)
    coef_y(6)=4*P1(2)-4*P4(2)+4*P5(2)-4*P6(2)

    coef_z(1)=P1(3)
    coef_z(2)=-3*P1(3)-P2(3)+4*P4(3)
    coef_z(3)=-3*P1(3)-P3(3)+4*P6(3)
    coef_z(4)=2*P1(3)+2*P2(3)-4*P4(3)
    coef_z(5)=2*P1(3)+2*P3(3)-4*P6(3)
    coef_z(6)=4*P1(3)-4*P4(3)+4*P5(3)-4*P6(3)

    do count=1,n_order
      F_x(count)=coef_x(1)+coef_x(2)*U(count)+coef_x(3)*V(count)+coef_x(4)&
          &*U(count)**2+coef_x(5)*V(count)**2+coef_x(6)*U(count)*V(count)
      F_y(count)=coef_y(1)+coef_y(2)*U(count)+coef_y(3)*V(count)+coef_y(4)&
          &*U(count)**2+coef_y(5)*V(count)**2+coef_y(6)*U(count)*V(count)
      F_z(count)=coef_z(1)+coef_z(2)*U(count)+coef_z(3)*V(count)+coef_z(4)&
          &*U(count)**2+coef_z(5)*V(count)**2+coef_z(6)*U(count)*V(count)
      U_x(count)=coef_x(2)+2*coef_x(4)*U(count)+coef_x(6)*V(count)
      U_y(count)=coef_y(2)+2*coef_y(4)*U(count)+coef_y(6)*V(count)
      U_z(count)=coef_z(2)+2*coef_z(4)*U(count)+coef_z(6)*V(count)
      V_x(count)=coef_x(3)+2*coef_x(5)*V(count)+coef_x(6)*U(count)
      V_y(count)=coef_y(3)+2*coef_y(5)*V(count)+coef_y(6)*U(count)
      V_z(count)=coef_z(3)+2*coef_z(5)*V(count)+coef_z(6)*U(count)
      nP_x(count)=U_y(count)*V_z(count)-U_z(count)*V_y(count);
      nP_y(count)=U_z(count)*V_x(count)-U_x(count)*V_z(count);
      nP_z(count)=U_x(count)*V_y(count)-U_y(count)*V_x(count);
      dS(count)=sqrt(nP_x(count)**2+nP_y(count)**2+nP_z(count)**2);
      nP_x(count)=nP_x(count)/dS(count)
      nP_y(count)=nP_y(count)/dS(count)
      nP_z(count)=nP_z(count)/dS(count)
    enddo
    return
  end subroutine eval_quadratic_patch_UV


*/



void cross_product_3d(double *x, double *y, double *z) {
  // this routine computes z = x X y, where X is the cross product

  z[0] = x[1]*y[2] - x[2]*y[3];
  z[1] = x[2]*y[0] - x[0]*y[3];
  z[2] = x[0]*y[1] - x[1]*y[0];

  return;
}


#include "cisurf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>





void read_msh( BaseMesh *meshout, long id, char *name, char *filename) {
  //
  // This routine reads in "filename" which is assumed to be a .msh
  // file and returns a data structure of mesh elements in meshout
  //
  FILE *fileptr;

  printf("\n\nLoading geometry file: %s\n", filename);
  fileptr = fopen( filename, "r+" );
  if (fileptr == NULL) {
    printf("The file is not opened properly, exiting.");
    exit(0);
  }

  // read in header info
  long aux1, aux2, aux3, nverts, nelems; 
  fscanf(fileptr, "%ld %ld %ld %ld %ld", &aux1, &aux2, &aux3, &nverts, &nelems);


  // assign some basic info to the mesh structure
  meshout->id = id;
  meshout->name = name;
  meshout->nverts = nverts;
  meshout->nelems = nelems;

  // allocate space for all the vertices and elements
  meshout->verts = (double *) malloc( 3*nverts*sizeof(double) );
  
  // load in the actual information
  long i, aux4, aux5, aux6, aux7, aux8;
  double x,y,z;
  for (i=0; i<nverts; i++) {
    fscanf(fileptr, "%ld %ld %ld %ld %ld %ld %ld %ld", &aux1, &aux2, &aux3,
           &aux4, &aux5, &aux6, &aux7, &aux8);
    fscanf(fileptr, "%lf %lf %lf", &x, &y, &z);

    meshout->verts[3*i] = x;
    meshout->verts[3*i+1] = y;
    meshout->verts[3*i+2] = z;
    
  }

  // now read in the actual elements
  fscanf(fileptr, "%ld", &aux1);
  long nv;

  meshout->elements = (BaseElement *) malloc(nelems*sizeof(BaseElement));
  nv = 6;
  for (i=0; i<nelems; i++) {
    fscanf(fileptr, "%ld %ld %ld %ld %ld %ld %ld %ld", &aux1, &aux2, &aux3,
           &aux4, &aux5, &aux6, &aux7, &aux8);
    fscanf(fileptr, "%ld %ld %ld %ld %ld %ld", &aux1, &aux2, &aux3,
           &aux4, &aux5, &aux6);

    meshout->elements[i].id = i;
    
    meshout->elements[i].nv = nv;
    meshout->elements[i].ivs = (long *) malloc( 6*sizeof(long) );
    meshout->elements[i].ivs[0] = aux1-1;
    meshout->elements[i].ivs[1] = aux2-1;
    meshout->elements[i].ivs[2] = aux3-1;
    meshout->elements[i].ivs[3] = aux4-1;
    meshout->elements[i].ivs[4] = aux5-1;
    meshout->elements[i].ivs[5] = aux6-1;

    meshout->elements[i].gtype = "tria2";
    
  }

  //printf("for example, gtype = %s\n", meshout->elements[50].gtype);
  //exit(0);


  // the data is loaded in, time to close the file
  fclose(fileptr);

  
  // compute the centroids and bounding spheres
  double *verts, *vert1, *vert2, *vert3, r, d;
  long j;
  verts = meshout->verts;
  
  for (i=0; i<nelems; i++) {
    vert1 = &verts[ meshout->elements[i].ivs[0] ];
    vert2 = &verts[ meshout->elements[i].ivs[1] ];
    vert3 = &verts[ meshout->elements[i].ivs[2] ];

    x = (vert1[0] + vert2[0] + vert3[0])/3;
    y = (vert1[1] + vert2[1] + vert3[1])/3;
    z = (vert1[2] + vert2[2] + vert3[2])/3;
    meshout->elements[i].centroid[0] = x;
    meshout->elements[i].centroid[1] = y;
    meshout->elements[i].centroid[2] = z;

    r = -1;
    d = sqrt( (vert1[0]-x)*(vert1[0]-x) + (vert1[1]-y)*(vert1[1]-y)
              + (vert1[2]-z)*(vert1[2]-z) );
    if (d > r) r = d;

    d = sqrt( (vert2[0]-x)*(vert2[0]-x) + (vert2[1]-y)*(vert2[1]-y)
              + (vert2[2]-z)*(vert2[2]-z) );
    if (d > r) r = d;

    d = sqrt( (vert3[0]-x)*(vert3[0]-x) + (vert3[1]-y)*(vert3[1]-y)
              + (vert3[2]-z)*(vert3[2]-z) );
    if (d > r) r = d;

    meshout->elements[i].radius = r;

  }

  //cisurf_print_element_info( &(meshout->elements[25]) );
  //exit(0);


  // compute the pseudo normal vectors used in constructing a smoothed mesh,
  // these are a function of the basemesh
  long *counter;
  counter = (long *) malloc( nverts*sizeof(long) );
  double *pns;

  meshout->pseudoNormals = (double *) malloc( nverts*3*sizeof(double) );
  pns = meshout->pseudoNormals;

  // initialize the counter and the pseudonormals vector
  for (i=0; i<nverts; i++) {
    counter[i] = 0;
    pns[3*i] = 0;
    pns[3*i+1] = 0;
    pns[3*i+2] = 0;
  }

  // now scan through each element, and for each vertex compute its normal
  // vector, add it to the pseduonormal vector and add one to the counter --
  // then at the end, divide by the counter to compute the averaged pseudonormal

  double uvs[100];
  uvs[0] = 0.0;
  uvs[1] = 0.0;

  uvs[2] = 1.0;
  uvs[3] = 0.0;

  uvs[4] = 0.0;
  uvs[5] = 1.0;


  uvs[0] = 0.5;
  uvs[1] = 0.0;

  uvs[2] = 0.5;
  uvs[3] = 0.5;

  uvs[4] = 0.0;
  uvs[5] = 0.5;


  long iv;
  double xyzs[1000], dus[1000], dvs[1000];
  double das[1000], normals[1000], verts2[1000];

  for (i=0; i<nelems; i++) {

    for (j=0; j<6; j++) {
      iv = meshout->elements[i].ivs[j];
      verts2[3*j] = verts[3*iv];
      verts2[3*j+1] = verts[3*iv+1];
      verts2[3*j+2] = verts[3*iv+2];
    }

    eval_tria2_pn( 6, uvs, verts2, xyzs, dus, dvs, das, normals );

    for (j=0; j<6; j++) {
      iv = meshout->elements[i].ivs[j];
      counter[iv] = counter[iv]+1;
      pns[3*iv] = pns[3*iv] + normals[3*j];
      pns[3*iv+1] = pns[3*iv+1] + normals[3*j+1];
      pns[3*iv+2] = pns[3*iv+2] + normals[3*j+2];
    }
  }

  // now scale by the counter
  for (i=0; i<nverts; i++) {
    pns[3*i] = pns[3*i]/counter[i];
    pns[3*i+1] = pns[3*i+1]/counter[i];
    pns[3*i+2] = pns[3*i+2]/counter[i];
  }


  free(counter);






  return;
}





void eval_tria2_pn( long n, double *uv, double *verts, double *xyz, double *du,
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
    cross_prod3d_(&du[3*i], &dv[3*i], &normal[3*i]);

    da[i] = sqrt( normal[3*i]*normal[3*i] + normal[3*i+1]*normal[3*i+1]
                  + normal[3*i+2]*normal[3*i+2] );
    normal[3*i] = normal[3*i]/(da[i]);
    normal[3*i+1] = normal[3*i+1]/(da[i]);
    normal[3*i+2] = normal[3*i+2]/(da[i]);

  }

  return;

}









void print_base_mesh_info( BaseMesh *mesh1, long iflong ) {
  //
  // print out the information contained in the mesh structure
  //
  printf("\n");
  printf("- - - - mesh information - - - - \n");
  printf("mesh.name   = %s\n", mesh1->name);
  printf("mesh.nverts = %ld\n", mesh1->nverts);
  printf("mesh.nelems = %ld\n", mesh1->nelems);

  if (iflong == 1) {

    printf("\n");

    long i;
    for (i=0; i<mesh1->nverts; i++) {
      printf("vertex %d   = (%e %e %e)\n", i, mesh1->verts[3*i],
             mesh1->verts[3*i+1], mesh1->verts[3*i+2] );
    }
  }

  printf("- - - end mesh information - - - \n");
  printf("\n");
  
}





void print_base_element_info( BaseElement *elem ) {
  //
  // print out the information contained in a single element
  //
  printf("\n");
  printf("- - - - element information - - - - \n");
  printf("element.id          = %ld\n", elem->id);
  printf("element.gtype       = %s\n", elem->gtype);
  printf("element.nv          = %ld\n", elem->nv);
  printf("element.ivs         =");

  long i;
  for (i=0; i<(elem->nv); i++) printf(" %ld", elem->ivs[i]);
  printf("\n");

  printf("element.centroid    = (%e, %e, %e)\n", elem->centroid[0],
         elem->centroid[1], elem->centroid[2]);
  printf("element.radius      = %e\n", elem->radius);
  printf("- - - end element information - - - \n");
  printf("\n");
  
}

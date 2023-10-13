
#include "cisurf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>





void readMSH(baseMesh *meshout, long id, char *name, char *filename) {
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

  meshout->elements = (baseElement *) malloc(nelems*sizeof(baseElement));
  nv = 6;
  for (i=0; i<nelems; i++) {
    fscanf(fileptr, "%ld %ld %ld %ld %ld %ld %ld %ld", &aux1, &aux2, &aux3,
	   &aux4, &aux5, &aux6, &aux7, &aux8);
    fscanf(fileptr, "%ld %ld %ld %ld %ld %ld", &aux1, &aux2, &aux3,
	   &aux4, &aux5, &aux6);

    meshout->elements[i].id = i;
    
    meshout->elements[i].nv = nv;
    meshout->elements[i].ivs = (long *) malloc( 6*sizeof(long) );
    meshout->elements[i].ivs[0] = aux1;
    meshout->elements[i].ivs[1] = aux2;
    meshout->elements[i].ivs[2] = aux3;
    meshout->elements[i].ivs[3] = aux4;
    meshout->elements[i].ivs[4] = aux5;
    meshout->elements[i].ivs[5] = aux6;

    meshout->elements[i].gtype = "tria2";
    
  }

  //printf("for example, gtype = %s\n", meshout->elements[50].gtype);
  //exit(0);

  
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
    
  fclose(fileptr);


  //cisurf_print_element_info( &(meshout->elements[25]) );
  //exit(0);
  
  return;
}





void printBaseMeshInfo( baseMesh *mesh1, long iflong ) {
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





void printBaseElementInfo( baseElement *elem ) {
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




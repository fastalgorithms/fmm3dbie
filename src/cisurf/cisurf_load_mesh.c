
#include "cisurf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>





void cisurf_read_msh(mesh *meshout, char *name, char *filename) {
  //
  // This routine reads in "filename" which is assumed to be a .msh
  // file and returns a data structure of mesh elements in meshout
  //
  FILE *fileptr;

  printf("Loading geometry file: %s\n", filename);
  fileptr = fopen( filename, "r+" );
  if (fileptr == NULL) {
    printf("The file is not opened properly, exiting.");
    exit(0);
  }

  // read in header info
  long aux1, aux2, aux3, nverts, nelems; 
  fscanf(fileptr, "%ld %ld %ld %ld %ld", &aux1, &aux2, &aux3, &nverts, &nelems);

  printf("after fscanf, aux1 = %ld\n", aux1);
  printf("after fscanf, aux2 = %ld\n", aux2);
  printf("after fscanf, aux3 = %ld\n", aux3);
  printf("after fscanf, nverts = %ld\n", nverts);
  printf("after fscanf, nelems = %ld\n", nelems);

  // assign some basic info to the mesh structure
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

  meshout->elements = (meshElement *) malloc(nelems*sizeof(meshElement));
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


  cisurf_print_element_info( &(meshout->elements[25]) );
  exit(0);
  
  return;
}





void cisurf_print_mesh_info( mesh *mesh1, long iflong ) {
  //
  // print out the information contained in the mesh structure
  //
  printf("\n");
  printf("- - - mesh information - - - \n");
  printf("mesh.name   = %s\n", mesh1->name);
  printf("mesh.nverts = %ld\n", mesh1->nverts);
  printf("mesh.nelems = %ld\n", mesh1->nelems);
  printf("\n");

  if (iflong == 1) {
  
    long i;
    for (i=0; i<mesh1->nverts; i++) {
      printf("vertex %d   = (%e %e %e)\n", i, mesh1->verts[3*i],
	     mesh1->verts[3*i+1], mesh1->verts[3*i+2] );
    }
  }
  
  printf("- - - end mesh information - - - \n");
  printf("\n");
  
}





void cisurf_print_element_info( meshElement *elem ) {
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



/* subroutine readmsh(Geometry1, filename, norder_skel, norder_smooth) */
/*   use ModType_Smooth_Surface */
/*   implicit none */

/*   !! This subroutine opens a msh file and loads the information in a */
/*   !! variable of type Geometry */

/*   ! */
/*   ! Input */
/*   !   filename - the file to read */
/*   !   norder_skel - order to discretize the skeleton patches */
/*   !   norder_smooth - order to discretize the smoothed patches */
/*   ! */
/*   ! Output */
/*   !   Geometry1 - data structure for geometry */
/*   ! */

/*   !List of calling arguments */
/*   type (Geometry), intent(inout) :: Geometry1 */
/*   character(len=100), intent(in) :: filename */
/*   integer :: n_order_sf, nsk, nsf */
/*   integer  :: norder_skel, norder_smooth */
  
/*   integer umio,i,m,N,j,aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8 */
/*   integer :: ierror */


/*   Geometry1%ifflat = 0 */

/*   open(UNIT=8, FILE=trim(filename), STATUS='OLD', ACTION='READ', IOSTAT=ierror) */
/*   read(8,*) aux1,aux2,aux3,m, N */

/*   print * */
/*   print * */
  
/*   write (6,*) 'loading file ', trim(filename) */
/*   write (13,*) 'loading file ', trim(filename) */

/*   call prinf('ntri = *', n, 1) */
/*   !call prinf('npoints = *', m, 1) */

  

/*   nsk = (norder_skel+1)*(norder_skel+2)/2 */
/*   !call prinf('nsk = *', nsk, 1) */
/*   !stop */

/*   call prinf('num points on skeleton mesh = *', nsk*n, 1) */
  
/*   Geometry1%norder_skel = norder_skel */
/*   Geometry1%nskel = nsk */

/*   !Geometry1%norder_smooth = norder_smooth */
/*   !Geometry1%nsmooth = n_order_sf */

/*   Geometry1%norder_smooth = norder_smooth */
/*   nsf = (norder_smooth+1)*(norder_smooth+2)/2 */
/*   Geometry1%nsmooth = nsf */

/*   call prinf('num points on smooth mesh = *', nsf*n, 1) */
  
/*   Geometry1%n_order_sf = nsf */
/*   Geometry1%npoints=m */
/*   Geometry1%ntri=N */
/*   Geometry1%ntri_sk=N */
/*   Geometry1%n_Sf_points=N*nsf */
/*   Geometry1%n_Sk_points=N*nsk */

/*   if (allocated(Geometry1%Points)) then */
/*     deallocate(Geometry1%Points) */
/*   endif */
/*   if (allocated(Geometry1%Tri)) then */
/*     deallocate(Geometry1%Tri) */
/*   endif */

/*   allocate(Geometry1%Points(3,m)) */
/*   allocate(Geometry1%Tri(6,N)) */

/*   do j=1,m */
/*     read(8,*) aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8 */
/*     read(8,*) Geometry1%Points(1,j),Geometry1%Points(2,j),Geometry1%Points(3,j) */
/*   enddo */

/*   read(8,*) aux1 */

/*   do j=1,N */
/*     read(8,*) aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8 */
/*     read(8,*) aux1,aux2,aux3,aux4,aux5,aux6 */
/*     Geometry1%Tri(1,j)=aux1 */
/*     Geometry1%Tri(2,j)=aux2 */
/*     Geometry1%Tri(3,j)=aux3 */
/*     Geometry1%Tri(4,j)=aux4 */
/*     Geometry1%Tri(5,j)=aux5 */
/*     Geometry1%Tri(6,j)=aux6 */
/*   enddo */
  
/*   close (8) */

/*   return */
/* end subroutine readmsh */


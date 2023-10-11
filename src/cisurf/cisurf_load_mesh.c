
#include "cisurf.h"
#include "cisurf_load_mesh.h"
#include <stdio.h>
#include <stdlib.h>

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

  // allocate space for all the vertices
  meshout->verts = (double *)malloc(3*nverts*sizeof(double));

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
  fclose(fileptr);


  
  return;
}





void cisurf_print_mesh_info( mesh *mesh1 ) {
  //
  // print out the information contained in the mesh structure
  //
  printf("\n");
  printf("- - - mesh information - - - \n");
  printf("mesh.name   = %s\n", mesh1->name);
  printf("mesh.nverts = %ld\n", mesh1->nverts);
  printf("mesh.nelems = %ld\n", mesh1->nelems);
  printf("\n");

  long i;
  for (i=0; i<mesh1->nverts; i++) {
    printf("vertex %d   = (%e %e %e)\n", i, mesh1->verts[3*i],
	   mesh1->verts[3*i+1], mesh1->verts[3*i+2] );
  }
  printf("- - - end mesh information - - - \n");
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


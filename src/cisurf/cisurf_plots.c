
#include "cisurf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>




void cisurf_plot_mesh_vtk( mesh *mesh1, char *filename ) {
  //
  // This routine creates an ASCII vtk file which can plot mesh1.
  //
  
  long nelems, nverts;
  nelems = mesh1->nelems;
  nverts = mesh1->nverts;
  
  meshElement *elem1;
  FILE *fptr;

  fptr = fopen( filename, "w" );
  fprintf(fptr, "# vtk DataFile Version 3.0\n");
  fprintf(fptr, "MESHNAME: %s\n", mesh1->name);
  fprintf(fptr, "ASCII\n");
  fprintf(fptr, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(fptr, "POINTS %ld float\n", nverts);

  long i;
  for (i=0; i<nverts; i++) {
    fprintf(fptr, "%e %e %e\n", mesh1->verts[3*i], mesh1->verts[3*i+1],
	    mesh1->verts[3*i+2]);
  }

  // compute total number of points needed across all elements,
  // allowing for various kinds of elements which have various numbers
  // of vertices defining them
  long ntot;
  ntot = 0;
  for (i=0; i<nelems; i++) {
    ntot = ntot + 1 + mesh1->elements[i].nv;
  }

  fprintf(fptr, "CELLS %ld %ld\n", nelems, ntot);

  for (i=0; i<nelems; i++) {


  }
  

  /* if(ifflat.eq.0) then */
  /*   do i=1,ntri */
  /*     write(iunit1,'(a,i9,i9,i9,i9,i9,i9)') "6 ", Geometry1%Tri(1,i)-1, & */
  /*      Geometry1%Tri(2,i)-1,Geometry1%Tri(3,i)-1, & */
  /*      Geometry1%Tri(4,i)-1,Geometry1%Tri(5,i)-1, Geometry1%Tri(6,i)-1 */
  /*   enddo */
  /*   write(iunit1,'(a,i9)') "CELL_TYPES ", ntri */
  /*   do ipatch = 1,ntri */
  /*     write(iunit1,'(a)') "22" */
  /*   end do */

  /* elseif(ifflat.eq.1) then */
  /*   do i=1,ntri */
  /*     write(iunit1,'(a,i9,i9,i9)') "3 ", Geometry1%Tri(1,i)-1, & */
  /*      Geometry1%Tri(2,i)-1,Geometry1%Tri(3,i)-1 */
  /*   enddo */
  /*   write(iunit1,'(a,i9)') "CELL_TYPES ", ntri */
  /*   do ipatch = 1,ntri */
  /*     write(iunit1,'(a)') "5" */
  /*   end do */
  /* endif */


  /* write(iunit1,'(a)') "" */
  /* write(iunit1,'(a,i9)') "POINT_DATA ", Geometry1%npoints */
  /* write(iunit1,'(a,i4)') "SCALARS normals float ", 3 */
  /* write(iunit1,'(a)') "LOOKUP_TABLE default" */
  /* do i = 1,Geometry1%npoints */
  /*   write(iunit1,'(E11.5,2x,E11.5,2x,e11.5)') & */
  /*     Geometry1%Normal_Vert(1,i),& */
  /*     Geometry1%Normal_Vert(2,i),& */
  /*     Geometry1%Normal_Vert(3,i) */
  /*   rr = Geometry1%Normal_Vert(1,i)**2 + & */
  /*        Geometry1%Normal_Vert(2,i)**2 + & */
  /*        Geometry1%Normal_Vert(3,i)**2 */
  /*   print *, i, rr */
  /* end do */

  /* close(iunit1) */


  fclose(fptr);

  

  return; 
}

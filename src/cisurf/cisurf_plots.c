
#include "cisurf.h"
#include "cprini_long.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>




void plot_base_mesh_vtk( BaseMesh *mesh1, char *filename, char *vertsfile ) {
  //
  // This routine creates an ASCII vtk file which can plot mesh1, and a separte
  // one that will plot the vertices used in describing the mesh.
  //


  long nelems, nverts;
  nelems = mesh1->nelems;
  nverts = mesh1->nverts;


  BaseElement *elem1;
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

  long j;
  for (i=0; i<nelems; i++) {
    fprintf(fptr, "%ld", mesh1->elements[i].nv);
    for (j=0; j<mesh1->elements[i].nv; j++) {
      fprintf(fptr, " %ld", mesh1->elements[i].ivs[j] );
    }
    fprintf(fptr, "\n");
  }

  // print the cell types
  fprintf(fptr, "CELL_TYPES %ld\n", nelems);
  for (i=0; i<nelems; i++) {
    fprintf(fptr, "22\n");
  }

  // plot the z value of each node so we can color the thing
  fprintf(fptr, "\n");
  fprintf(fptr, "POINT_DATA %ld\n", nverts);
  fprintf(fptr, "SCALARS z-value float 1\n");
  fprintf(fptr, "LOOKUP_TABLE default\n");

  for (i=0; i<nverts; i++) {
    fprintf(fptr, "%e\n", mesh1->verts[3*i+2]);
  }

  fclose(fptr);




  // and now plot the vertices along with pseudonormal information
  fptr = fopen( vertsfile, "w" );
  fprintf(fptr, "# vtk DataFile Version 3.0\n");
  fprintf(fptr, "VERTSNAME: %s\n", mesh1->name);
  fprintf(fptr, "ASCII\n");
  fprintf(fptr, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(fptr, "POINTS %ld float\n", nverts);

  for (i=0; i<nverts; i++) {
    fprintf(fptr, "%e %e %e\n", mesh1->verts[3*i], mesh1->verts[3*i+1],
            mesh1->verts[3*i+2]);
  }

  fprintf(fptr, "\n" );
  fprintf(fptr, "CELLS %ld %ld\n", nverts, 2*nverts);

  for (i=0; i<nverts; i++) {
    fprintf(fptr, "%ld %ld\n", 1, i);
  }

  // print the cell types
  fprintf(fptr, "\n" );
  fprintf(fptr, "CELL_TYPES %ld\n", nverts);
  for (i=0; i<nverts; i++) {
    fprintf(fptr, "1\n");
  }

  // plot the z value of each node so we can color the thing
  fprintf(fptr, "\n");
  fprintf(fptr, "POINT_DATA %ld\n", nverts);
  fprintf(fptr, "SCALARS z-value float 1\n");
  fprintf(fptr, "LOOKUP_TABLE default\n");

  for (i=0; i<nverts; i++) {
    fprintf(fptr, "%e\n", mesh1->verts[3*i+2]);
  }

  fprintf(fptr, "\n");
  fprintf(fptr, "VECTORS pseudoNormals float\n");

  for (i=0; i<nverts; i++) {
    fprintf(fptr, "%e %e %e\n", mesh1->pseudoNormals[3*i],
            mesh1->pseudoNormals[3*i+1], mesh1->pseudoNormals[3*i+2] );
    //fprintf(fptr, "0.0 0.0 1.0\n" );
  }

  fclose(fptr);



  return;
}






void plot_skeleton_mesh_vtk( SkelMesh *mesh1, char *meshfile, char *nodefile ) {
  //
  // This routine creates an ASCII vtk file which can plot mesh1, which is of
  // type SkelMesh and contains quadrature nodes, normals, weights, etc. The
  // element as well as the quadrature nodes are plotted.
  //
  // Input:
  //   mesh1 - a skeleton mesh to plot
  //   meshfile - the file to dump the mesh info, e.g. "skelmesh.vtk"
  //   nodefile - a separate file to dump the discretization nodes along with
  //     other info like du, dv, normals, etc., e.g. "skelnodes.vtk"
  //
  // The mesh and the node info are exported in two separate files, which must
  // be loaded separately.
  //

  long nelems;
  nelems = mesh1->nelems;

  SkelElement *elem1;
  FILE *fptr;

  // dump out the points and various scalars and values

  fptr = fopen( nodefile, "w" );
  fprintf(fptr, "# vtk DataFile Version 3.0\n");
  fprintf(fptr, "MESHNAME: %s\n", mesh1->name);
  fprintf(fptr, "ASCII\n");
  fprintf(fptr, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(fptr, "\n");

  long npts = 0;
  long i,j, npols;
  //for (i=0; i<nelems; i++) {
  //  npts = npts + mesh1->elements[i].npols;
  // }

  npts = 3*nelems;

  cprinf("inside plot skelmesh, nelems = ", &nelems, 1);

  fprintf(fptr, "POINTS %ld float\n", npts);
  for (i=0; i<nelems; i++) {
    //npols = mesh1->elements[i].npols;
    //for (j=0; j<npols; j++) {
    //fprintf(fptr, "%e %e %e\n", mesh1->elements[i].srcvals[j].xyz[0],
    //          mesh1->elements[i].srcvals[j].xyz[1],
    //          mesh1->elements[i].srcvals[j].xyz[2] );
    fprintf(fptr, "%e %e %e\n", mesh1->elements[i].verts[0],
            mesh1->elements[i].verts[1],
            mesh1->elements[i].verts[2] );
    fprintf(fptr, "%e %e %e\n", mesh1->elements[i].verts[3],
            mesh1->elements[i].verts[4],
            mesh1->elements[i].verts[5] );
    fprintf(fptr, "%e %e %e\n", mesh1->elements[i].verts[6],
            mesh1->elements[i].verts[7],
            mesh1->elements[i].verts[8] );
      //}
  }



  /* long i; */
  /* for (i=0; i<nverts; i++) { */
  /*   fprintf(fptr, "%e %e %e\n", mesh1->verts[3*i], mesh1->verts[3*i+1], */
  /*           mesh1->verts[3*i+2]); */
  /* } */

  /* // compute total number of points needed across all elements, */
  /* // allowing for various kinds of elements which have various numbers */
  /* // of vertices defining them */
  /* long ntot; */
  /* ntot = 0; */
  /* for (i=0; i<nelems; i++) { */
  /*   ntot = ntot + 1 + mesh1->elements[i].nv; */
  /* } */

  /* fprintf(fptr, "CELLS %ld %ld\n", nelems, ntot); */

  /* long j; */
  /* for (i=0; i<nelems; i++) { */
  /*   fprintf(fptr, "%ld", mesh1->elements[i].nv); */
  /*   for (j=0; j<mesh1->elements[i].nv; j++) { */
  /*     fprintf(fptr, " %ld", mesh1->elements[i].ivs[j]-1); */
  /*   } */
  /*   fprintf(fptr, "\n"); */
  /* } */

  /* // print the cell types */
  /* fprintf(fptr, "CELL_TYPES %ld\n", nelems); */
  /* for (i=0; i<nelems; i++) { */
  /*   fprintf(fptr, "22\n"); */
  /* } */

  /* // plot the z value of each node so we can color the thing */
  /* fprintf(fptr, "\n"); */
  /* fprintf(fptr, "POINT_DATA %ld\n", nverts); */
  /* fprintf(fptr, "SCALARS z-value float 1\n"); */
  /* fprintf(fptr, "LOOKUP_TABLE default\n"); */

  /* for (i=0; i<nverts; i++) { */
  /*   fprintf(fptr, "%e\n", mesh1->verts[3*i+2]); */
  /* } */

  /* /\* elseif(ifflat.eq.1) then *\/ */
  /* /\*   do i=1,ntri *\/ */
  /* /\*     write(iunit1,'(a,i9,i9,i9)') "3 ", Geometry1%Tri(1,i)-1, & *\/ */
  /* /\*      Geometry1%Tri(2,i)-1,Geometry1%Tri(3,i)-1 *\/ */
  /* /\*   enddo *\/ */
  /* /\*   write(iunit1,'(a,i9)') "CELL_TYPES ", ntri *\/ */
  /* /\*   do ipatch = 1,ntri *\/ */
  /* /\*     write(iunit1,'(a)') "5" *\/ */
  /* /\*   end do *\/ */
  /* /\* endif *\/ */


  fclose(fptr);

  return;
}

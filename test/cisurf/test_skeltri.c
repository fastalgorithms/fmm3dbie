
//
// test program for loading various kinds of mesh files, and making a
// plot of them
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cisurf.h"
#include "cprini_long.h"

int main(int argc, char **argv) {



  //printf("number of arguments received = %d\n", argc);

  //long i;
  //for (i=0; i<argc; i++) {
  //  printf("argument %d is %s\n", i+1, argv[i]);
  //}

  // load in a mesh file and construct a skeleton on it

  cprin_init("cfort.13", "stdout");

  char *filename;
  
  if (argc == 1) {
    filename = "../../geometries/msh/round_2nd_tri.msh";
  }
  else {
    filename = argv[1];
  }
  
  //printf("\n\n");
  cprin_skipline(2);
  cprin_message( "File to be loaded: " );
  cprin_message( filename );
  cprin_skipline(2);
  //printf("File to be loaded: %s\n", filename);

  
  BaseMesh basemesh1;
  long id;
  id = 1;
  read_msh( &basemesh1, id, "base mesh", filename);

  // print the mesh infor
  print_base_mesh_info( &basemesh1, 0);

  // dump the mesh as a vtk file
  plot_base_mesh_vtk( &basemesh1, "basemesh.vtk", "baseverts.vtk" );

  exit(0);

  // now construct a skeleton mesh
  SkelMesh skelmesh1;
  skelmesh1.id = 0;
  skelmesh1.name = "skeleton mesh";

  long norder;
  norder = 8;
  create_skeleton( &basemesh1, &skelmesh1, norder );

  plot_skeleton_mesh_vtk( &skelmesh1, "skelmesh.vtk", "skelnodes.vtk" );

  return 0;
}


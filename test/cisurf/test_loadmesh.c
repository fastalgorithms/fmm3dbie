
//
// test program for loading various kinds of mesh files, and making a
// plot of them
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cisurf.h"

int main(int argc, char **argv) {



  //printf("number of arguments received = %d\n", argc);

  //long i;
  //for (i=0; i<argc; i++) {
  //  printf("argument %d is %s\n", i+1, argv[i]);
  //}

  char *filename;
  
  if (argc == 1) {
    filename = "../../geometries/msh/round_2nd_tri.msh";
  }
  else {
    filename = argv[1];
  }
  
  printf("\n\n");
  printf("File to be loaded: %s\n", filename);

  
  baseMesh mesh1;
  long id;
  id = 1;
  readMSH(&mesh1, id, "base mesh", filename);

  // print the mesh infor
  printBaseMeshInfo(&mesh1, 0);

  // dump the mesh as a vtk file
  plotBaseMeshVTK( &mesh1, "basemesh.vtk" );
  
  return 0;
}



//
// test program for loading various kinds of mesh files, and making a
// plot of them
//

#include <stdio.h>
#include <string.h>
#include "cisurf.h"

int main() {


  char filename[] = "../../geometries/msh/round_2nd_tri.msh";

  printf("File to be loaded: %s\n", filename);

  mesh mesh1;
  cisurf_read_msh(&mesh1, "a rounded avocado thing", filename);

  // print the mesh infor
  cisurf_print_mesh_info(&mesh1, 1);

  // dump the mesh as a vtk file
  cisurf_plot_mesh_vtk( &mesh1, "basemesh.vtk" );
  
  return 0;
}

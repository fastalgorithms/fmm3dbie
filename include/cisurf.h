
//
// header file for cisurf routines, containing data structures for the
// mesh, skeleton mesh, etc.
//
// (c) Mike O'Neil 2023
//     moneil@flatironinstitute.org
//


typedef struct cisurfMeshElement {

  // an index number of the element
  long id;
  
  // the element type, given as a string of maximum length 8, supported types are:
  //   "tria1" - flat triangles
  //   "tria2" - quadrature triangles, specified by vertices and midpoints
  //   "quad1" - flat quadrilateral
  //   "quad2" - quadratic quadrilateral, specified by vertices and midpoints
  char *gtype;

  // a collection of the vertices needed to define the element,
  // oriented counterclockwise around the element so as to give the
  // proper outward normal -- this varies, and will depend on the
  // element type
  long nv;
  long *ivs;
  
  
  // some basic measurement info of the element
  double centroid[3];
  double radius;

} meshElement;






typedef struct cisurfMesh {

  // name the mesh?
  char *name;
  
  // total number of vertices
  long nverts;
  double *verts;
  
  // total number of elements, and an array of them
  long nelems;
  meshElement *elements;
  
} mesh;





// function prototypes

void cisurf_read_msh(mesh *meshout, char *name, char *filename);

void cisurf_print_mesh_info( mesh *mesh1, long iflong);

void cisurf_print_element_info( meshElement *elem );

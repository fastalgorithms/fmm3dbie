
//
// header file for cisurf routines, containing data structures for the
// bash mesh, skeleton mesh, surface mesh, etc., and various function
// prototypes
//
// (c) Mike O'Neil 2023
//     moneil@flatironinstitute.org
//


typedef struct baseElement {

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

} baseElement;






typedef struct baseMesh {

  // name the mesh?
  long id;
  char *name;
  
  // total number of vertices
  long nverts;
  double *verts;
  
  // total number of elements, and an array of them
  long nelems;
  baseElement *elements;
  
} baseMesh;






// function prototypes

// - - - mesh loading routines - - -
void readMSH(baseMesh *meshout, char *name, char *filename);

void printBaseMeshInfo( baseMesh *mesh1, long iflong);

void printBaseElementInfo( baseElement *elem );



// - - - plotting routines - - -
void plotBaseMeshVTK( baseMesh *mesh1, char *filename );

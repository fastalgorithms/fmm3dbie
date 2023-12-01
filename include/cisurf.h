
//
// header file for cisurf routines, containing data structures for the
// bash mesh, skeleton mesh, surface mesh, etc., and various function
// prototypes
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





typedef struct pointInfo {
  // stores xyz coordinates, surface differentials, and normal
  // information at *some* point. It is assumed that the normal is
  // actually normalized to have unit length, but obviously du and dv
  // don't.
  double xyz[3];
  double du[3];
  double dv[3];
  double normal[3];
  
} pointInfo;





typedef struct coefsInfo {
  // stores the coefficient expansions for a skeleton element, may
  // correspond to Legendre or Koornwinder, the data structure is the
  // same
  double *x, *y, *z;
  double *dxdu, *dydu, *dzdu;
  double *dxdv, *dydv, *dzdv;

} coefsInfo;




typedef struct skelElement {

  // the struct storing information for a single element on the
  // skeleton mesh, this means quadrature nodes, etc.
  long id;

  // gytype can have different values:
  //   "tria1" - flat triangles
  //   "tria2" - quadrature triangles, specified by vertices and midpoints
  //   "quad1" - flat quadrilateral
  //   "quad2" - quadratic quadrilateral, specified by vertices and midpoints
  char *gtype;

  // collection of vertices needed to define the element, same as for
  // a base element
  long nv, *ivs;

  // discretization info
  long norder;

  pointInfo *srcvals;
  coefsInfo *coefs;

  // some basic measurement info
  double centroid[3];
  double radius;

  
} skelElement;





typedef struct skelMesh {

  // a collection of skeleton elements, much in the style of base
  // meshes but with additional info
  long id;
  char *name;

  // total number of vertices
  long nverts;
  double *verts;

  // total number of elements, and an array of them
  long nelems;
  skelElement *elements;


} skelMesh;





// function prototypes

// - - - mesh loading routines - - -
void readMSH(baseMesh *meshout, long id, char *name, char *filename);

void printBaseMeshInfo( baseMesh *mesh1, long iflong);

void printBaseElementInfo( baseElement *elem );


// - - - skeleton construction routines - - -
void makeSkeleton( baseMesh *baseMesh1, skelMesh *skelMesh1, long norder );

// - - - plotting routines - - -
void plotBaseMeshVTK( baseMesh *mesh1, char *filename );



// external functions, maybe in fortran
void get_vioreanu_nodes_wts_( int *norder, int *npols, double *uvs,
                              double *whts );

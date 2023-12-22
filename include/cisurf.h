
//
// header file for cisurf routines, containing data structures for the
// bash mesh, skeleton mesh, surface mesh, etc., and various function
// prototypes
//


typedef struct BaseElement {

  // an index number of the element
  long id;
  
  // the element type, given as a string of maximum length 8, supported types are:
  //   "tria1" - flat triangles
  //   "tria2" - quadratic triangles, specified by vertices and midpoints
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

} BaseElement;






typedef struct BaseMesh {

  // name the mesh?
  long id;
  char *name;
  
  // vertex information
  long nverts;
  double *verts;
  double *pseudoNormals;
  
  // total number of elements, and an array of them
  long nelems;
  BaseElement *elements;
  
} BaseMesh;





typedef struct PointInfo {
  // stores xyz coordinates, surface differentials, and normal information at
  // *some* point on the patch. It is assumed that the normal is actually
  // normalized to have unit length, but obviously du and dv don't.
  double xyz[3];
  double du[3];
  double dv[3];
  double normal[3];
  double pseudoNormal[3];

} PointInfo;





typedef struct CoefsInfo {
  // stores the coefficient expansions for a skeleton element, may
  // correspond to Legendre or Koornwinder, the data structure is the
  // same
  double *x, *y, *z;
  double *dxdu, *dydu, *dzdu;
  double *dxdv, *dydv, *dzdv;

} CoefsInfo;




typedef struct SkelElement {

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
  long nv;
  double *verts;

  // discretization info
  //   norder - the order at which the element is discretized
  //   npols - the number of discretization nodes
  long norder;
  long npols;
  double *whts;

  PointInfo *srcvals;
  CoefsInfo *coefs;

  // some basic measurement info
  double centroid[3];
  double radius;

  
} SkelElement;





typedef struct SkelMesh {

  // a collection of skeleton elements, much in the style of base
  // meshes but with additional info
  long id;
  char *name;

  // total number of elements, and an array of them
  long nelems;
  SkelElement *elements;


} SkelMesh;





// function prototypes

// - - - mesh loading routines - - -
void read_msh(BaseMesh *meshout, long id, char *name, char *filename);

void print_base_mesh_info( BaseMesh *mesh1, long iflong);

void print_base_element_info( BaseElement *elem );




// - - - skeleton mesh routines - - -
void create_skeleton( BaseMesh *basemesh1, SkelMesh *skelmesh1, long norder );

void eval_tria2( long n, double *uv, double *verts, double *xyz, double *du,
                 double *dv, double *da, double *normal );

void eval_tria2_pn( long n, double *uv, double *verts, double *xyz, double *du,
                 double *dv, double *da, double *normal );

void cross_product_3d(double *x, double *y, double *z);

void print_skeleton_element_info( SkelElement *element );


// - - - plotting routines - - -
void plot_base_mesh_vtk( BaseMesh *mesh1, char *meshfiles, char *vertsfile );

void plot_skeleton_mesh_vtk( SkelMesh *mesh1, char *meshfile, char *nodefile );

// external functions, maybe in fortran
void get_vioreanu_nodes_wts_( long *norder, long *npols, double *uvs,
                              double *whts );

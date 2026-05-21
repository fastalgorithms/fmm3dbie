#include <stdio.h>
#include <string.h>
#include <stddef.h>

#ifndef MWF77_RETURN
#define MWF77_RETURN int
#endif

#if defined(MWF77_CAPS)
#define MWF77_multiscale_mesher MULTISCALE_MESHER_UNIF_REFINE
#elif defined(MWF77_UNDERSCORE1)
#define MWF77_multiscale_mesher multiscale_mesher_unif_refine_ 
#elif defined(MWF77_UNDERSCORE0)
#define MWF77_multiscale_mesher multiscale_mesher_unif_refine 
#else
#define MWF77_multiscale_mesher multiscale_mesher_unif_refine__
#endif

void multiscale_mesher_unif_refine_cfname_(char *, int64_t*, int64_t*, char *, int64_t*, int64_t *, int64_t *, int64_t *, double *, char *, int64_t *);
void f2cstr_(char *);

#ifdef __cplusplus
extern "C"
#endif

int main(int argc, char **argv)
{
  char *filenamein;
  filenamein = "../../geometries/meshes/sphere.msh";

  char *filenameout;
  filenameout = "../../geometries/sphere";

  int64_t ifiletype = 1;
  int64_t norder_skel = 12;
  int64_t norder_smooth = 4;
  int64_t nrefine = 1;
  int64_t adapt_flag = 1;
  double rlam = 2.5;
  int64_t ier = 0;
  int64_t ifcad = 0;
  char *fcad;
  fcad = "tmp";
  

  multiscale_mesher_unif_refine_cfname_(filenamein, &ifiletype, &ifcad, fcad, &norder_skel, &norder_smooth, 
     &nrefine, &adapt_flag, &rlam, filenameout, &ier);
  return 0;

}

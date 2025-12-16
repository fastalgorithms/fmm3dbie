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

void multiscale_mesher_unif_refine_cfname_(char *, int*, int*, int *, int *, int *, double *, char *, int *);
void f2cstr_(char *);

#ifdef __cplusplus
extern "C"
#endif

int main(int argc, char **argv)
{
  char *filenamein;
  filenamein = "/Users/mrachh/git/fmm3dbie/geometries/meshes/sphere.msh";

  char *filenameout;
  filenameout = "/Users/mrachh/git/fmm3dbie/geometries/sphere";

  int ifiletype = 1;
  int norder_skel = 12;
  int norder_smooth = 4;
  int nrefine = 1;
  int adapt_flag = 1;
  double rlam = 2.5;
  int ier = 0;
  

  multiscale_mesher_unif_refine_cfname_(filenamein, &ifiletype, &norder_skel, &norder_smooth, 
     &nrefine, &adapt_flag, &rlam, filenameout, &ier);
  return 0;

}


#include <string.h>
#include "cprini_long.h"

//
//******************************************************************************
//
// This is the end of the debugging code, and the beginning of the printing
// functions proper. The code below contains five user-callable functions. Their
// functions are as follows.
//
//  cprin_init - communicates to the function cprin8 (NOT a user-callable
//        function) the names of the two files on which the output will be
//        written. Please note that the name "stdout" is treated differently
//        from all other names: specifying it will cause the output to be
//        written on the standard output - screen, most probably.
//
//  cprin_start_stop - starts and stops printing to either or both output files
//
//  cprin_message - prints a text message
//
//  cprin_skipline - skips a certain number of lines in the printout
//
//  cprinz - prints out complex numbers
//
//  cprinz_matrix - prints out a complex valued matrix
//
//  cprinf - prints a long integer array
//
//  cprind - print out double array
//
//  cprind_matrix - print out a double valued matrix
//
//******************************************************************************
//


void cprin_init(char *str17,  char *str27) {
  char *mes;
  float *ap;
  long *afp;
  double *adp;
  char *acp;
  long n, m;
  long itype;
  long i1, i2;
  //long *adp,*ap,itype,i1,i2,mes,afp,n,acp;

  /* This function initializes which files to print to.  */
  /* To print to ONLY one file, set one of str17, str18 to " ".  */
  /* i.e., make the call cprini("stdout"," ") to print ONLY to  */
  /* the screen.  Anything other than " ", such as "    ", will cause   */
  /* a file to be created, in this case "    ".  */

  itype=11;
  cprin_master(mes, ap, afp, adp, acp, m, n, itype, str17, str27, i1, i2);
  return;
}





void cprin_start_stop(long i1, long i2)
{
  char *mes;
  float *ap;
  long *afp;
  double *adp;
  char *acp;
  long n, m;
  long itype;
  char *str17,*str27;

  // . . . stop/resume.  i1 and i2 control printing to file
  // str17 and str27, respectively.  0 means stop printing,
  // and 1 means print.

  itype=21;
  cprin_master(mes,ap,afp,adp,acp,m,n,itype,str17,str27,i1,i2);
  return;
}







void cprin_message(char *mes)
{
  float *ap;
  long *afp;
  char *acp;
  long itype;
  char *str17, *str27;
  long i1, i2;
  double *adp;
  long m, n;
  
  itype = 99;
  cprin_master(mes,ap,afp,adp,acp,m,n,itype,str17,str27,i1,i2);
  return;
}





void cprin_skipline(long n)
{
  float *ap;
  long *afp;
  char *acp;
  long itype;
  char *str17, *str27;
  long i1, i2;
  long status;
  char *mes = NULL;
  double *adp;
  long m;
  
  //  printf("in cprin_skipline, n  = %d\n", n);
  //exit(status);
  
  //long *ap,*afp,itype,str17,str27,i1,i2,acp;
  /*
   *   Print double precision data
   */
  itype = 77;
  
  //printf("in cprin_skipline, itype  = %d\n", itype);
  
  cprin_master(mes,ap,afp,adp,acp,m,n,itype,str17,str27,i1,i2);
  return;
}





void cprinz(char *mes, double _Complex *adp, long n)
{
  float *ap;
  long *afp;
  char *acp;
  long itype;
  char *str17, *str27;
  long i1, i2;
  long status;
  long n2, m;

  //printf("n  = %ld\n", n);
  //exit(status);
  
  //long *ap,*afp,itype,str17,str27,i1,i2,acp;
  /*
   *   Print double precision data
   */
  itype=7;
  n2 = 2*n;
  cprin_master(mes,ap,afp, (double *)adp,
               acp,m,n2,itype,str17,str27,i1,i2);
  return;
}





void cprinz_matrix(char *mes, double _Complex *adp, long m, long n)
{
  // print out a double complex valued matrix with m rows and n columns, the
  // assumption is that the data is stored in column major format in adp

  float *ap;
  long *afp;
  char *acp;
  long itype;
  char *str17, *str27;
  long i1, i2;
  long status;
  long n2;

  itype=33;
  n2 = 2*n;
  cprin_master(mes,ap,afp, (double *)adp,
               acp,m,n2,itype,str17,str27,i1,i2);
  return;
}





void cprinf(char *mes, long *afp, long n)
{
  float *ap;
  long itype, i1, i2, m;
  double *adp;
  char *acp, *str17, *str27;

  //
  // print long integer data
  //
  itype = 2;
  cprin_master(mes, ap, afp, adp, acp, m, n, itype, str17, str27, i1, i2);
  return;
}





void cprind(char *mes, double *adp, long n) {

  float *ap;
  long itype, i1, i2, m, *afp;
  char *acp, *str17, *str27;
  
  //
  // print double precision data
  //
  itype = 3;
  cprin_master(mes, ap, afp, adp, acp, m, n, itype, str17, str27, i1, i2);
  return;
}





void cprind_matrix(char *mes, double *adp, long m, long n)
{
  float *ap;
  long *afp;
  char *acp;
  long itype;
  char *str17, *str27;
  long i1, i2;
  long status;
  long n2;

  //printf("n  = %ld\n", n);
  //exit(status);
  
  //long *ap,*afp,itype,str17,str27,i1,i2,acp;
  /*
   *   Print double precision data
   */
  itype=33;
  cprin_master(mes,ap,afp,adp,acp,m,n,itype,str17,str27,i1,i2);
  return;
}





void cprin_master(char *mes, float *ap, long *afp, double *adp, char *acp, long m,
                  long n, long itype, char *str17, char *str27, long i1, long i2)
{
  static long ifprint1,ifprint2,ifstr1,ifstr2;
  static FILE *st1,*st2;
  long iii;
  
  // If this is the initialization call - open the 
  // files str17, str27 after checking if either is null, i.e. " "

  //printf("in cprin_master, itype = %d\n", itype);
  //return;

  
  if(itype==11)
    {
      ifprint1=0;
      ifprint2=0;
      
      ifstr1=0;
      ifstr2=0;
      
      iii=strcmp(str17," ");
      if(iii != 0) ifprint1=1;
      if(iii != 0) ifstr1=1;
      
      iii=strcmp(str27," ");
      if(iii != 0) ifprint2=1;
      if(iii != 0) ifstr2=1;
      
      if(ifprint1 == 1)
        {
          iii=strcmp(str17,"stdout");
          if(iii != 0) st1=fopen(str17,"w");
          if(iii == 0) st1=stdout;
        }
      
      if(ifprint2 == 1) {
        iii=strcmp(str27,"stdout");
        if(iii != 0) st2=fopen(str27,"w");
        if(iii == 0) st2=stdout;
      }
      
      return;
    }
  
  //
  //   If this is the "stop/resume" call - stop/resume printing 
  //
  if(itype==21){
    if(i1==0) ifprint1=0;
    if((i1 != 0) && (ifstr1 != 0)) ifprint1=1;
    if(i2==0) ifprint2=0;
    if((i2 != 0) && (ifstr2 != 0)) ifprint2=1;
    return;
  }

  //printf("ifprint1 = %d\n", ifprint1);
  //printf("ifprint2 = %d\n", ifprint2);

  
  if(ifprint1 !=0) cprin_all(mes,ap,afp,adp,acp,m,n,itype,st1);
  if(ifprint2 !=0) cprin_all(mes,ap,afp,adp,acp,m,n,itype,st2);
  
  return;
}





void cprin_all(char *mes, float *ap, long *afp, double *adp, char *acp,
               long m, long n, long itype, FILE *str)
{
  long i;
  //
  // Process the message
  //

  if(mes) fprintf(str,"%s\n",mes);
  if (itype == 99) return;


  
  
  // Process the double precision data to be printed
  if(itype == 3) {
    for(i=0; i<n; i=i+1) {
      fprintf(str,"  %11.4le", adp[i]);
      if(i%6==5 || i==n-1) fprintf(str,"\n");
    }
    return;
  }


  //
  // print out the double precision matrix
  //
  if(itype ==33) {
    long ijk = 0;
    long j;
    
    //printf("m = %d\n", m);
    //printf("n = %d\n", n);      
    //return;

    
    for(i=0; i<m; i=i+1) {
      for (j=0; j<n; j++) {
        fprintf(str,"  %11.4le", adp[ijk]);
        ijk++;
      }
      fprintf(str,"\n");
    }
    return;
  }

  


  /*
   *   Process the complex double precision data to be printed
   */
  if(itype ==7)
    {
      for(i=0; i<n; i=i+2)
        {
          fprintf(str," (%11.4le %+11.4le)", adp[i], adp[i+1]);
          if(i%6==4 || i==n-2) fprintf(str,"\n");
        }
      return;
    }

/*
*   Process the integer data to be printed
*/
    if(itype ==2)
    {
        for(i=0; i<n; i=i+1)
        {
            fprintf(str," %7d", afp[i]);
            if(i%10==9 || i==n-1) fprintf(str,"\n");
        }
        return;
    }

/*
*   Process the single precision data to be printed
*/
    if(itype ==1)
    {
        for(i=0; i<n; i=i+1)
        {
            fprintf(str,"  %11.4le", ap[i]);
            if(i%6==5 || i==n-1) fprintf(str,"\n");
        }
        return;
    }

/*
*   Process the character data to be printed
*/
    if(itype==4) {
      for(i=0; i<n; i++) {
        fprintf(str,"%c", acp[i]);
        if(i%60==59 || i==n-1) fprintf(str,"\n");
      }
      return;
    }

    
    //
    // insert a line skip
    //
    if(itype == 77) {
      for (i=0; i<n; i++) {
        fprintf(str,"\n");
      }
    }

    return;
}

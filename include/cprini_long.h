#ifndef CPRINI_LONG_H_
#define CPRINI_LONG_H_


#include <stdio.h>


void cprin_init(char *str1, char *str2);

void cprin_master(char *mes, float *ap, long *afp, double *adp, char *acp,
                  long m, long n, long itype, char *str17, char *str27,
                  long i1, long i2);

void cprin_all(char *mes, float *ap, long *afp, double *adp, char *acp,
               long m, long n, long itype, FILE *str);

void cprinf(char *mes, long *ip, long n);

void cprind(char *mes, double *adp, long n);

void cprind_matrix(char *mes, double *adp, long m, long n);

void cprinz(char *mes, double _Complex *adp, long n);

void cprin_message(char *mes);

void cprin_skipline(long n);

void cprin_start_stop(long i1, long i2);



#endif // CPRINI_LONG_H_

#pragma once
#include "bvh_binary.h"

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void get_near_pairs_(const int32_t* ns, const double* sources, const int32_t* nt, const double* targets,
                     int32_t* n_pairs, int32_t* near_pairs);

void free_near_pairs_(int32_t* n_pairs, int32_t* near_pairs);

#ifdef __cplusplus
}
#endif

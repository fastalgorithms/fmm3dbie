#include "bvh_fortran_interface.h"

#ifdef __cplusplus
extern "C" {
#endif

void get_near_pairs_(const int32_t* ns, const double* sources, const int32_t* nt, const double* targets,
                     int32_t* n_pairs, int32_t* near_pairs) {
    // set lbvh boxes and build bvh
    LBVH lbvh;
#pragma omp parallel for
    for (int32_t i = 0; i < ns[0]; ++i) {
        int32_t idx = 6 * i;
        lbvh.pos.push_back(glm::vec3((float)sources[idx + 0], (float)sources[idx + 1], (float)sources[idx + 2]));
        lbvh.pos.push_back(glm::vec3((float)sources[idx + 3], (float)sources[idx + 4], (float)sources[idx + 5]));
    }
    lbvh.build();

    // set overlap query
    std::vector<AABB> aabbs;
#pragma omp parallel for
    for (int32_t i = 0; i < nt[0]; ++i) {
        AABB aabb;
        int32_t idx = 6 * i;
        aabb.min = glm::vec3((float)targets[idx + 0], (float)targets[idx + 1], (float)targets[idx + 2]);
        aabb.max = glm::vec3((float)targets[idx + 3], (float)targets[idx + 4], (float)targets[idx + 5]);
    }
    std::vector<std::pair<uint32_t, uint32_t>> col_list;

    // compute overlap
    lbvh.overlap_aabbs(aabbs, col_list);

    // return near pairs
    n_pairs[0] = (int32_t)col_list.size();
    near_pairs = (int32_t*)malloc(n_pairs[0] * 2 * sizeof(int32_t));

#pragma omp parallel for
    for (int32_t i = 0; i < n_pairs[0]; ++i) {
        near_pairs[2 * i + 0] = col_list[i].first;
        near_pairs[2 * i + 1] = col_list[i].second;
    }
}

void free_near_pairs_(int32_t* n_pairs, int32_t* near_pairs) {
    n_pairs[0] = 0;
    free(near_pairs);
}

#ifdef __cplusplus
}
#endif

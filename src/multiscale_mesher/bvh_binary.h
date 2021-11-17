#pragma once

#include "glm/glm.hpp"

#include <algorithm> // std::min, std::sort
#include <atomic>    // std::atomic
#include <cstring>
#include <omp.h>
#include <vector> // std::vector

const auto maximum = std::numeric_limits<float>::max();
const auto minimum = std::numeric_limits<float>::lowest();
const auto invalid = std::numeric_limits<uint32_t>::max();

inline static uint32_t to_uint(float f) {
    uint32_t ui;
    memcpy(&ui, &f, sizeof(float));
    return ui;
}

inline static float to_float(uint32_t ui) {
    float f;
    memcpy(&f, &ui, sizeof(uint32_t));
    return f;
}

struct AABB {
    union {
        std::atomic<uint32_t> atomic_min[3];
        glm::vec3 min;
    };
    union {
        std::atomic<uint32_t> atomic_max[3];
        glm::vec3 max;
    };

    inline void update_min(std::atomic<uint32_t>& min_value, const uint32_t value) noexcept {
        uint32_t old_value = min_value;
        while (old_value > value && !min_value.compare_exchange_weak(old_value, value)) {
        }
    }

    inline void update_max(std::atomic<uint32_t>& max_value, const uint32_t value) noexcept {
        uint32_t old_value = max_value;
        while (old_value < value && !max_value.compare_exchange_weak(old_value, value)) {
        }
    }

    inline void update_min(glm::vec3 val) noexcept { min = glm::min(val, min); }

    inline void update_max(glm::vec3 val) noexcept { max = glm::max(val, max); }

    AABB() {}

    void init() {
        min = glm::vec3(maximum);
        max = glm::vec3(minimum);
    }

    AABB(const AABB& box) {
        min = box.min;
        max = box.max;
    }

    void operator=(const AABB& box) {
        min = box.min;
        max = box.max;
    }

    // expansion by point doesn't need atomic ops
    inline void expand(const glm::vec3& p) {
        min = glm::min(p, min);
        max = glm::max(p, max);
    }

    inline void expand(const AABB& box) {
        for (auto i = 0; i < 3; ++i) {
            update_min(atomic_min[i], box.atomic_min[i]);
            update_max(atomic_max[i], box.atomic_max[i]);
        }
    }

    inline glm::vec3 normalize(const glm::vec3& p) const { return (p - min) / (max - min); }

    inline float halved_surf_area() const {
        const auto w = glm::max(glm::vec3(0), max - min);
        return (w[0] * w[1] + w[1] * w[2] + w[2] * w[0]);
    }

    inline bool contain_point(const glm::vec3& p) const {
        return glm::all(glm::lessThanEqual(min, p)) && glm::all(glm::greaterThanEqual(max, p));
    }

    inline bool overlap_aabb(const AABB& aabb) const {
        for (int d = 0; d < 3; ++d)
            if (min[d] > aabb.max[d] || max[d] < aabb.min[d])
                return false;
        return true;
    }
};

struct alignas(8) Node {
    AABB box;

    // Child Indices
    // if lowest bit is 1, leaf
    //                  0, node

    uint32_t l_child;
    uint32_t r_child;

    Node() : l_child(invalid), r_child(invalid) {
        box.atomic_max[0] = 0;
        box.atomic_max[1] = 0;
        box.atomic_max[2] = 0;

        box.atomic_min[0] = to_uint(maximum);
        box.atomic_min[1] = to_uint(maximum);
        box.atomic_min[2] = to_uint(maximum);
    }

    inline bool is_l_leaf() const noexcept {
        const unsigned char mask = 1;
        return (l_child & mask);
    }

    inline bool is_r_leaf() const noexcept {
        const unsigned char mask = 1;
        return (r_child & mask);
    }
};

// Ciprian Apetrei "Fast and Simple Agglomerative LBVH Construction"
// http://diglib.eg.org/handle/10.2312/cgvc.20141206.041-044
struct LBVH {
    // Thinking Parallel, Part III: Tree Construction on the GPU
    // https://devblogs.nvidia.com/thinking-parallel-part-iii-tree-construction-gpu/

    // Expands a 10-bit integer into 30 bits
    // by inserting 2 zeros after each bit.
    inline uint32_t expand_bits(uint32_t v) {
        v = (v * 0x00010001u) & 0xFF0000FFu;
        v = (v * 0x00000101u) & 0x0F00F00Fu;
        v = (v * 0x00000011u) & 0xC30C30C3u;
        v = (v * 0x00000005u) & 0x49249249u;
        return v;
    }

    // Calculates a 30-bit Morton code for the
    // given 3D point located within the unit cube [0,1].
    inline uint32_t morton_3d(float x, float y, float z) {
        x = std::min(std::max(x * 1024.0f, 0.0f), 1023.0f);
        y = std::min(std::max(y * 1024.0f, 0.0f), 1023.0f);
        z = std::min(std::max(z * 1024.0f, 0.0f), 1023.0f);
        uint32_t xx = expand_bits((uint32_t)x);
        uint32_t yy = expand_bits((uint32_t)y);
        uint32_t zz = expand_bits((uint32_t)z);
        return xx * 4 + yy * 2 + zz;
    }

    inline bool is_leaf(const uint32_t node_id) const noexcept {
        const unsigned char mask = 1;
        return (node_id & mask);
    }

    AABB box;      // Scene Bound
    uint32_t root; // Root Node ID

    std::vector<Node> nodes;    // LBVH nodes (num_boxes - 1)
    std::vector<glm::vec3> pos; // min max of bboxes

    void build();

    // col_list stores the leaf nodes overlapping with aabb
    void overlap_aabb(const AABB& aabb, std::vector<uint32_t>& col_list);
    void overlap_aabbs(const std::vector<AABB>& aabbs, std::vector<std::pair<uint32_t, uint32_t>>& col_list);
    bool overlap_lchild(const AABB& aabb, const uint32_t& node);
    bool overlap_rchild(const AABB& aabb, const uint32_t& node);
    inline bool is_l_leaf(const uint32_t& node);
    inline bool is_r_leaf(const uint32_t& node);
    inline uint32_t l_index(const uint32_t& node);
    inline uint32_t r_index(const uint32_t& node);
    double nodes_halved_surf_area();
};

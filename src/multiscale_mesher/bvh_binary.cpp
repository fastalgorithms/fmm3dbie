#include "bvh_binary.h"

#include <chrono>    // std::chrono
#include <execution> // std::execution::par
#include <iomanip>   // std::setprecision
#include <iostream>  // std::cout

// delta function in sec3 of the paper
// "Fast and Simple Agglomerative LBVH Construction"
inline uint32_t delta(const std::vector<glm::uvec3>& leaves, const uint32_t id) {
    return leaves[id + 1].z ^ leaves[id].z;
}

void LBVH::build() {
    const auto start = std::chrono::steady_clock::now();

    // number of boxes
    const auto num_boxes = pos.size() / 2;

    // object bounding box
    box.init();
    for (uint32_t i = 0; i < num_boxes; ++i) {
        box.update_min(pos[2 * i + 0]);
        box.update_max(pos[2 * i + 1]);
    }

    // allocate pair <reference, morton code>
    std::vector<glm::uvec3> leaves(num_boxes);

    // set order
#pragma omp parallel for
    for (uint32_t i = 0; i < num_boxes; ++i)
        leaves[i].x = i;

    // set morton code
    std::for_each(std::execution::par, leaves.begin(), leaves.end(), [&](glm::uvec3& leaf) {
        const auto id0 = leaf.x * 2 + 0;
        const auto id1 = leaf.x * 2 + 1;
        const auto centroid = (pos[id0] + pos[id1]) / 2.0f;
        const auto unitcube = box.normalize(centroid);          // coord in unit cube
        leaf.z = morton_3d(unitcube.x, unitcube.y, unitcube.z); // morton code
        // std::cout << "morton: " << leaf.z << ", x: " << unitcube.x << ", y: " << unitcube.y
        //          << ", z: " << unitcube.z << "\n";
    });

    // sort leaves by morton code in ascending order
    std::sort(std::execution::par, std::begin(leaves), std::end(leaves),
              [&](const glm::uvec3& l, const glm::uvec3& r) { return (l.z < r.z); });

#pragma omp parallel for
    for (uint32_t i = 0; i < num_boxes; ++i)
        leaves[i].y = i;

    // number of nodes
    const auto num_nodes = num_boxes - 1;

    // allocate inner nodes
    nodes.resize(num_nodes);

    // otherBounds in algorithm 1 of the paper
    // "Massively Parallel Construction of Radix Tree Forests for the Efficient Sampling of Discrete
    // Probability Distributions" https://arxiv.org/pdf/1901.05423.pdf
    std::vector<std::atomic<uint32_t>> other_bounds(num_nodes);
    std::for_each(std::execution::par, other_bounds.begin(), other_bounds.end(),
                  [&](std::atomic<uint32_t>& b) { b.store(invalid); });

    // for each leaf
    std::for_each(std::execution::par, leaves.begin(), leaves.end(), [&](glm::uvec3& leaf) {
        // current leaf/node id
        auto current = leaf.y;

        // range
        uint32_t left = current;
        uint32_t right = current;

        // leaf aabb
        const auto id = leaf.x * 2;
        AABB aabb;
        aabb.init();
        aabb.min = pos[id + 0];
        aabb.max = pos[id + 1];

        // current is leaf or node?
        auto is_leaf = true;
        while (1) {
            // the whole range is covered
            if (left == 0 && right == num_nodes) {
                root = current;
                break;
            }

            // leaf/node index
            const auto index = is_leaf ? leaves[current].x * 2 + 1 : current * 2;

            // choose parent
            uint32_t previous, parent;
            if (0 == left || (right != num_nodes && delta(leaves, right) < delta(leaves, left - 1))) {
                // parent is right and "left" doesn't change
                parent = right;
                previous = other_bounds[parent].exchange(left);
                if (invalid != previous)
                    right = previous;
                nodes[parent].l_child = index;
            } else {
                // parent is left and "right" doesn't change
                parent = left - 1;
                previous = other_bounds[parent].exchange(right);
                if (invalid != previous)
                    left = previous;
                nodes[parent].r_child = index;
            }

            // expand aabb
            nodes[parent].box.expand(aabb);

            // terminate this thread
            if (invalid == previous)
                break;

            assert(left < right);

            // ascend
            current = parent;
            aabb = nodes[current].box;
            is_leaf = false;
        }
    });

//    std::for_each(std::execution::par, nodes.begin(), nodes.end(), [&](Node& node) {
//        node.box.min += box.min;
//        node.box.max += box.min;
//    });

    const auto end = std::chrono::steady_clock::now();
    std::cout << "bvh: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"
              << std::endl;

    // debug
    double sah = nodes_halved_surf_area();
    std::cout << "number nodes " << nodes.size() << std::endl;
    std::cout << "number leaves " << leaves.size() << std::endl;
    std::cout << "sah: " << std::setprecision(16) << sah << std::endl;
}

void LBVH::overlap_aabb(const AABB& aabb, std::vector<uint32_t>& col_list) {
    uint32_t stack[64];
    uint32_t* stack_ptr = stack;
    *stack_ptr++ = invalid;

    // buggy, node must be internal node for now, to fix
    uint32_t node = root;
    do {
        bool overlap_l = overlap_lchild(aabb, node);
        bool overlap_r = overlap_rchild(aabb, node);

        bool is_l_leaf = this->is_l_leaf(node);
        bool is_r_leaf = this->is_r_leaf(node);

        if (overlap_l && is_l_leaf)
            col_list.push_back(l_index(node));
        if (overlap_r && is_r_leaf)
            col_list.push_back(r_index(node));

        bool traverse_l = (overlap_l && !is_l_leaf);
        bool traverse_r = (overlap_r && !is_r_leaf);

        if (!traverse_l && !traverse_r) {
            node = *--stack_ptr;
        } else {
            node = (traverse_l) ? (nodes[node].l_child >> 1) : (nodes[node].r_child >> 1);
            if (traverse_l && traverse_r)
                *stack_ptr++ = (nodes[node].r_child >> 1);
        }
    } while (node != invalid);
}

void LBVH::overlap_aabbs(const std::vector<AABB>& aabbs, std::vector<std::pair<uint32_t, uint32_t>>& col_list) {
    // work space for each aabb
    std::vector<uint32_t> list_entry;
    list_entry.reserve(10);
    std::vector<std::vector<uint32_t>> aabb_col_list(aabbs.size(), list_entry);

    std::for_each(std::execution::par, aabbs.begin(), aabbs.end(), [&](const AABB& aabb) {
        uint32_t idx = (uint32_t)(&aabb - &aabbs[0]);
        overlap_aabb(aabb, aabb_col_list[idx]);
    });

    for (uint32_t i = 0; i < aabb_col_list.size(); ++i) {
        for (const auto& col : aabb_col_list[i]) {
            col_list.push_back(std::make_pair(i, col));
        }
    }
}

bool LBVH::overlap_lchild(const AABB& aabb, const uint32_t& node) {
    uint32_t lchild = l_index(node);
    if (is_l_leaf(node)) {
        AABB leaf_aabb;
        leaf_aabb.init();
        leaf_aabb.min = pos[2 * lchild + 0];
        leaf_aabb.max = pos[2 * lchild + 1];
        return aabb.overlap_aabb(leaf_aabb);
    } else {
        return aabb.overlap_aabb(nodes[lchild].box);
    }
}

bool LBVH::overlap_rchild(const AABB& aabb, const uint32_t& node) {
    uint32_t rchild = r_index(node);
    if (is_r_leaf(node)) {
        AABB leaf_aabb;
        leaf_aabb.init();
        leaf_aabb.min = pos[2 * rchild + 0];
        leaf_aabb.max = pos[2 * rchild + 1];
        return aabb.overlap_aabb(leaf_aabb);
    } else {
        return aabb.overlap_aabb(nodes[rchild].box);
    }
}

inline bool LBVH::is_l_leaf(const uint32_t& node) { return nodes[node].is_l_leaf(); }

inline bool LBVH::is_r_leaf(const uint32_t& node) { return nodes[node].is_r_leaf(); }

inline uint32_t LBVH::l_index(const uint32_t& node) { return (nodes[node].l_child >> 1); }

inline uint32_t LBVH::r_index(const uint32_t& node) { return (nodes[node].r_child >> 1); }

double LBVH::nodes_halved_surf_area() {
    double sah = 0.0;
    for (const auto& n : nodes) {
        sah += n.box.halved_surf_area();
    }
    return sah;
}

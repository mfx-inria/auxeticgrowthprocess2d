#ifndef CONNECTIVITY_H
#define CONNECTIVITY_H

#include <vector>
#include <queue>
#include <string>
#include <limits>

const unsigned int solid_phase = 0;
const unsigned int void_phase = 1;

inline int mod(int x, int m) {
    int r = x % m;
    return r < 0 ? r + m : r;
}

unsigned int num_connected_components_image_phase(const std::vector<unsigned int>& image, const unsigned int dim_x, const unsigned int dim_y, bool is_periodic, unsigned int phase) {
    /**
     Given a labeling, returns the number of phase connected components (4-neighbourhood, given phase)
     Optionally in a periodic domain.
     **/
    unsigned int num_connected_components = 0;
    std::vector<unsigned int> image_visited(dim_x * dim_y);
    std::fill(image_visited.begin(), image_visited.end(), 0);
    for (unsigned int x = 0; x < dim_x; x++) {
        for (unsigned int y = 0; y < dim_y; y++) {
            unsigned int p = x + y * dim_x;
            if (image_visited[p] == 0) {
                if (image[p] == phase) {
                    num_connected_components++;
                    std::queue<unsigned int> to_visit;
                    to_visit.push(p);
                    while (!to_visit.empty()) {
                        unsigned int pv = to_visit.front();
                        unsigned int xv = pv % dim_x;
                        unsigned int yv = pv / dim_x;
                        to_visit.pop();
                        if (image_visited[pv] == 0) {
                            image_visited[pv] = 1;
                            for (int i = -1; i <= 1; i = i + 1) {
                                for (int j = -1; j <= 1; j = j + 1) {
                                    if (abs(i) + abs(j) == 1) {
                                        int xn = xv + i;
                                        if (is_periodic) {
                                            xn = mod(xn, dim_x);
                                        }
                                        if (xn >= 0 && xn < static_cast<int>(dim_x)) {
                                            int yn = yv + j;
                                            if (is_periodic) {
                                                yn = mod(yn, dim_y);
                                            }
                                            if (yn >= 0 && yn < static_cast<int>(dim_y)) {
                                                unsigned int pn = xn + yn * dim_x;
                                                if (image[pn] == phase) {
                                                    to_visit.push(pn);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return num_connected_components;
}

#endif

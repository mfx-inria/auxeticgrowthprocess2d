#ifndef GROWTH_H
#define GROWTH_H

#include <set>
#include <queue>
#include <vector>
#include <limits>
#include <string>

const int neigh_indices[] = {0, 1, 0, -1, 1, 0, -1, 0};  // 4-connected neighbourhood indices

struct GrowthInfo {
    double distance;
    unsigned int s;
    int x, y, osx, osy;
};

void init_growth_info(GrowthInfo& gi, double distance, unsigned int s, int x, int y, int osx, int osy) {
    gi.distance = distance;
    gi.s = s; gi.x = x; gi.y = y; gi.osx = osx; gi.osy = osy;
}

class CompareGrowthInfo {
 public:
    bool operator()(GrowthInfo n1, GrowthInfo n2) const {
        if (n1.distance > n2.distance) {
            return true;
        } else {
            if (n1.distance < n2.distance) return false;
            return n1.s < n2.s;  // lexicographic comparison when at equal distance
        }
    }
};

inline void funcmod(int& a, int& o, const int p) {
    while (a >= p) {
        a = a - p;
        o -= 1;
    }
    while (a < 0) {
        a = a + p;
        o += 1;
    }
}

struct growth_stats_struct {
    unsigned int num_iterations_growth;
    double avg_size_queue, avg_iterations_loop_discard;
};

inline void update_cumulative_average(double& cma, unsigned int& n, const double val) {
    cma = (val + n * cma) / static_cast<double>(n + 1);
    n++;
}

void compute_discrete_growth(unsigned int image_size, std::vector<Vec2d>& sites, std::vector<unsigned int>& image_closest_site, const double max_growth_length_const, starshaped_struct& starshaped, starshaped_struct& max_growth_length_starshaped, const std::string& name, unsigned int num_growth_plots, std::vector<std::vector<GrowthInfo> >& Cw, growth_stats_struct& growth_stats) {
    std::vector<int> image_closest_site_osx(image_closest_site.size());  // auxiliary vectors
    std::vector<int> image_closest_site_osy(image_closest_site.size());
    std::priority_queue<GrowthInfo, std::vector<GrowthInfo >, CompareGrowthInfo> queue;  // priority queue
    // nucleation
    for (unsigned int s = 0; s < sites.size(); s++) {
        const Vec2d grid_pos = sites[s] * static_cast<double>(image_size);
        const unsigned int x = static_cast<unsigned int>(std::floor(grid_pos.x));
        const unsigned int y = static_cast<unsigned int>(std::floor(grid_pos.y));
        Vec2d coordinate = get_coordinate(x, y, image_size);
        GrowthInfo queue_element;
        init_growth_info(queue_element, get_distance_squared(coordinate, sites[s], starshaped), s, x, y, 0, 0);
        queue.push(queue_element);
    }
    // growth
    unsigned int num_closest_sites_labeled = 0, current_growth_plot = 0;
    growth_stats.num_iterations_growth = 0;
    growth_stats.avg_size_queue = growth_stats.avg_iterations_loop_discard = 0;
    while (!queue.empty()) {
        // extract and remove queue element with the lowest distance
        GrowthInfo queue_element = queue.top();
        queue.pop();
        update_cumulative_average(growth_stats.avg_size_queue, growth_stats.num_iterations_growth, static_cast<double>(queue.size()));
        // check if this element will be included in the cell, or will be discarded
        const unsigned int pos_element = queue_element.x + queue_element.y * image_size;
        if (image_closest_site[pos_element] == std::numeric_limits<unsigned int>::max()) {
            double max_growth_length = max_growth_length_const;
            if (!max_growth_length_starshaped.radial_spans.empty()) {  // variable max growth length case
                Vec2d site = sites[queue_element.s] + Vec2d(static_cast<double>(queue_element.osx), static_cast<double>(queue_element.osy));
                Vec2d direction = get_coordinate(queue_element.x, queue_element.y, image_size) - site;
                max_growth_length = get_interpolated_radial_span(direction, max_growth_length_starshaped);
            }
            const int discrete_max_growth_length = std::round(max_growth_length);
            const double max_growth_length_squared = max_growth_length * max_growth_length;
            const int image_size_int = static_cast<int>(image_size);
            bool discard = false;
            unsigned int iterations_loop_discard = 0;
            for (int xo = -discrete_max_growth_length; (xo <= discrete_max_growth_length) && !discard; ++xo) {
                int xp = queue_element.x + xo;
                int osx = queue_element.osx;
                funcmod(xp, osx, image_size_int);
                for (int yo = -discrete_max_growth_length; (yo <= discrete_max_growth_length) && !discard; ++yo) {
                    const Vec2d offset_length(static_cast<double>(xo), static_cast<double>(yo));
                    iterations_loop_discard++;
                    if (offset_length.length_squared() <= max_growth_length_squared) {
                        int yp = queue_element.y + yo;
                        int osy = queue_element.osy;
                        funcmod(yp, osy, image_size_int);
                        const unsigned int p = xp + yp * image_size;
                        if (image_closest_site[p] != std::numeric_limits<unsigned int>::max()) {
                            if (image_closest_site[p] != queue_element.s) {
                                discard = true;
                            } else if (image_closest_site_osx[p] != osx) {
                                discard = true;
                            } else if (image_closest_site_osy[p] != osy) {
                                discard = true;
                            }
                        }
                    }
                }
            }
            update_cumulative_average(growth_stats.avg_iterations_loop_discard, growth_stats.num_iterations_growth, iterations_loop_discard);
            if (!discard) {
                image_closest_site[pos_element] = queue_element.s;
                image_closest_site_osx[pos_element] = queue_element.osx;
                image_closest_site_osy[pos_element] = queue_element.osy;
                Cw[queue_element.s].push_back(queue_element);
                if (num_growth_plots > 0) {
                    const double relative_growth = static_cast<double>(num_closest_sites_labeled) / static_cast<double>(image_size * image_size);
                    if (static_cast<unsigned int>(std::floor(relative_growth * static_cast<double>(num_growth_plots))) >= current_growth_plot) {
                        if (current_growth_plot > 0) {
                            std::cout << "- Plotting porous material at relative growth = " << relative_growth * 100.0 << "%" << std::endl;
                            std::string str_relative_growth_percent = std::to_string(static_cast<unsigned int>(relative_growth * 100.0));
                            std::string filename_porous_material_png = name + "_" + str_relative_growth_percent + ".png";
                            plot_porous_material_and_sites_png(sites, image_size, image_size, filename_porous_material_png, image_closest_site, true, true);
                        }
                        current_growth_plot++;
                    }
                }
                // insert future candidate points of growth around a discrete local neighborhood (4-connected)
                for (unsigned int k = 0; k < 4; k++) {
                    int xn = queue_element.x + neigh_indices[k * 2];
                    int osx = queue_element.osx;
                    funcmod(xn, osx, image_size_int);
                    int yn = queue_element.y + neigh_indices[k * 2 + 1];
                    int osy = queue_element.osy;
                    funcmod(yn, osy, image_size_int);
                    if (image_closest_site[xn + yn * image_size] == std::numeric_limits<unsigned int>::max()) {
                        Vec2d coordinate = get_coordinate(xn, yn, image_size);
                        Vec2d site = sites[queue_element.s] + Vec2d(static_cast<double>(osx), static_cast<double>(osy));
                        GrowthInfo queue_element_neigh;
                        init_growth_info(queue_element_neigh, get_distance_squared(coordinate, site, starshaped), queue_element.s, xn, yn, osx, osy);
                        queue.push(queue_element_neigh);
                    }
                }
                num_closest_sites_labeled++;
            }
        }
        growth_stats.num_iterations_growth++;
    }
}

#endif

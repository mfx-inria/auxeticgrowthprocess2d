#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>

#include "vec2.h"
#include "starshaped.h"
#include "connectivity.h"
#include "plots.h"
#include "growth.h"
#include "point_process.h"
#include "parameters.h"
#include "regularization.h"


void initialize_starshaped(starshaped_struct& starshaped, SymmetryType& symmetry_type, unsigned int symmetry_degree, const std::string& filename_plot, bool plot_ppm, bool plot_pdf, bool is_starshaped_distance, bool plot_symmetry_axes, bool plot_interpolation_points, bool plot_interpolated_starshaped_set, double input_max_radial_span) {
     // initialize star-shaped set
    starshaped.symmetry_degree = symmetry_degree;
    if (symmetry_type == NoSymmetry) {
        starshaped.symmetry_degree = 0;
    }
    symmetrize_radial_spans(starshaped, symmetry_type, symmetry_degree);
    update_min_max_radial_spans(starshaped);
    starshaped.max_radial_span = std::max(starshaped.max_radial_span, input_max_radial_span);
    const unsigned int plot_size = 1000;
    if (plot_ppm) {
        std::string filename_plot_ppm = filename_plot + ".ppm";
        plot_starshaped_distance(plot_size, starshaped, filename_plot_ppm);
    }
    if (plot_pdf) {
        std::string filename_plot_pdf = filename_plot + ".pdf";
        plot_starshaped_distance_cairo(plot_size, starshaped, filename_plot_pdf, is_starshaped_distance, plot_symmetry_axes, plot_interpolation_points, plot_interpolated_starshaped_set);
    }
}


int main(int argc, char **argv) {
    // read input parameters
    if (argc != 2) {
        std::cout << "* Usage ./growthprocess2d [parameters file]" << std::endl;
        return -1;
    }
    parameters_struct parameters;
    if (!initialize_parameters(std::string(argv[1]), parameters)) return -1;

    std::string filename_distance_plot = parameters.name + "_distance";
    initialize_starshaped(parameters.starshaped, parameters.symmetry_type, parameters.symmetry_degree, filename_distance_plot, parameters.plot_starshaped_ppm, parameters.plot_starshaped_pdf, true, parameters.plot_symmetry_axes, parameters.plot_interpolation_points, parameters.plot_interpolated_starshaped_set, parameters.max_radial_span);
    if (!parameters.max_growth_length_starshaped.radial_spans.empty()) {
        std::string filename_max_growth_plot = parameters.name + "_max_growth";
        initialize_starshaped(parameters.max_growth_length_starshaped, parameters.symmetry_type, parameters.symmetry_degree, filename_max_growth_plot, parameters.plot_starshaped_ppm, parameters.plot_starshaped_pdf, false, parameters.plot_symmetry_axes, parameters.plot_interpolation_points, parameters.plot_interpolated_starshaped_set, parameters.max_length_growth_radial_span);
    }

    std::cout << "*** Point process..." << std::endl;
    std::vector<Vec2d> sites;
    if (parameters.point_process == Random) {
        double area_square = static_cast<double>(parameters.image_size * parameters.image_size);
        poisson_point_process(
            sites, parameters.num_sites_per_pixel, area_square, parameters.random_seed);
    } else if (parameters.point_process == RSA) {
        const unsigned int max_trials_rsa = 10000000;
        double rsa_max_dist = parameters.rsa_max_dist_rel_pixel / static_cast<double>(parameters.image_size);
        rsa_point_process(sites, rsa_max_dist, true, parameters.random_seed, max_trials_rsa);
        std::cout << "* Approximate RSA coverage = " << sites.size() * M_PI * rsa_max_dist * rsa_max_dist * 0.25 << std::endl;
    } else if (parameters.point_process == File) {
        file_point_process(sites, parameters.filename_points);
    }
    unsigned int num_sites = sites.size();
    std::cout << "* num_sites = " << num_sites << std::endl;
    if (num_sites == 0) {
        std::cerr << "[error] no sites" << std::endl;
        return -1;
    }

    std::vector<unsigned int> image_closest_site(parameters.image_size * parameters.image_size);
    std::fill(image_closest_site.begin(), image_closest_site.end(), std::numeric_limits<unsigned int>::max());

    std::cout << "*** Growth process..." << std::endl;
    auto start_growth_process = std::chrono::steady_clock::now();
    std::vector<std::vector<GrowthInfo>> Cw(sites.size());
    growth_stats_struct growth_stats;
    compute_discrete_growth(parameters.image_size, sites, image_closest_site, parameters.max_growth_length, parameters.starshaped, parameters.max_growth_length_starshaped, parameters.name, parameters.num_growth_plots, Cw, growth_stats);
    auto end_growth_process = std::chrono::steady_clock::now();
    unsigned int time_growth_process_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_growth_process - start_growth_process).count();
    std::cout << "* Elapsed time growth process: " << time_growth_process_ms << " ms" << std::endl;

    // compute statistics growth
    unsigned int dim_x = parameters.image_size;
    unsigned int dim_y = parameters.image_size;
    std::vector<unsigned int> image_porous_material(dim_x * dim_y);
    unsigned int num_void_pixels = 0;
    for (unsigned int y = 0; y < dim_y; y++) {
        for (unsigned int x = 0; x < dim_x ; x++) {
            unsigned int cu = void_phase;
            if (image_closest_site[x + y * dim_x] == std::numeric_limits<unsigned int>::max()) {
                cu = solid_phase;
            } else {
                num_void_pixels++;
            }
            image_porous_material[x + y * dim_x] = cu;
        }
    }

    double porosity = static_cast<double>(num_void_pixels) / static_cast<double>(dim_x * dim_y);
    std::cout << "* porosity = " << porosity << std::endl;
    unsigned int num_connected_components = num_connected_components_image_phase(image_porous_material, dim_x, dim_y, true, solid_phase);
    std::cout << "* num_connected_components=" << num_connected_components << std::endl;

    std::cout << "*** Cell regularization..." << std::endl;
    auto start_cell_regularization = std::chrono::steady_clock::now();
    unsigned int num_void_pixels_regularization = 0;
    std::vector<unsigned int> image_porous_material_regularized(dim_x * dim_y);
    regularize_cells(Cw, image_porous_material_regularized, dim_x, dim_y, num_void_pixels_regularization);
    auto end_cell_regularization = std::chrono::steady_clock::now();
    unsigned int time_cell_regularization_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_cell_regularization - start_cell_regularization).count();
    std::cout << "* Elapsed time cell regularization: " << time_cell_regularization_ms << " ms" << std::endl;

    double porosity_regularized = static_cast<double>(num_void_pixels_regularization) / static_cast<double>(dim_x * dim_y);
    std::cout << "* porosity_regularized = " << porosity_regularized << std::endl;
    unsigned int num_connected_components_regularized = num_connected_components_image_phase(image_porous_material_regularized, dim_x, dim_y, true, solid_phase);
    std::cout << "* num_connected_components_regularized=" << num_connected_components_regularized << std::endl;

    if (parameters.save_results_txt) {
        std::ofstream results_file((parameters.name + ".res").c_str());
        results_file << "num_sites=" << num_sites << std::endl;
        results_file << "porosity=" << porosity << std::endl;
        results_file << "porosity_regularized=" << porosity_regularized << std::endl;
        results_file << "num_connected_components=" << num_connected_components << std::endl;
        results_file << "num_connected_components_regularized=" << num_connected_components_regularized << std::endl;
        results_file << "time_growth_process_ms=" << time_growth_process_ms << std::endl;
        results_file << "time_cell_regularization_ms=" << time_cell_regularization_ms << std::endl;
        results_file << "num_iterations_growth=" << growth_stats.num_iterations_growth << std::endl;
        results_file << "avg_size_queue=" << growth_stats.avg_size_queue << std::endl;
        results_file << "avg_iterations_loop_discard=" << growth_stats.avg_iterations_loop_discard << std::endl;
        results_file.close();
    }
    if (parameters.save_porous_material_ppm) {
        std::string filename_porous_material_ppm = parameters.name + ".ppm";
        save_porous_material_ppm(filename_porous_material_ppm, image_porous_material_regularized, dim_x, dim_y);
    }
    if (parameters.plot_sites_pdf) {
        std::string filename_plot_sites =  parameters.name + "_point_process.pdf";
        plot_sites_pdf(sites, dim_x, dim_y, filename_plot_sites);
    }
    if (parameters.plot_porous_material_sites_png) {
        std::string filename_porous_material_png =  parameters.name + "_porous_material.png";
        plot_porous_material_and_sites_png(sites, dim_x, dim_y, filename_porous_material_png, image_porous_material_regularized, false, parameters.plot_sites_png);
        if (parameters.save_non_regularized) {
            std::string filename_porous_material_non_regularized_png =  parameters.name + "_porous_material_non_regularized.png";
            plot_porous_material_and_sites_png(sites, dim_x, dim_y, filename_porous_material_non_regularized_png, image_porous_material, false, parameters.plot_sites_png);
        }
    }
}

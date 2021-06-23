#ifndef PLOTS_H
#define PLOTS_H

#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>

#include <string>
#include <vector>
#include <limits>

const double color_point_process[3] = {136.0 / 255.0, 54.0 / 255.0, 54.0 / 255.0};
const double color_porous_material_growing[3] = {0.6, 0.6, 0.6};
const double color_cell_growing[3] = {1.0, 1.0, 1.0};
const double color_porous_material[3] = {0.0, 0.0, 0.0};
const double color_starshaped_distance[3] = {0.4, 0.6, 0.8};
const double color_euclidean_neighbourhood[3] = {0.4, 0.8, 0.6};
const double starshaped_distance_cairo_length_factor = 0.98;


void plot_starshaped_distance(unsigned int plot_size, starshaped_struct& starshaped, std::string& filename) {
    std::ofstream ofs;
    ofs.open(filename.c_str(), std::ios::binary);
    Vec2d origin(0.0, 0.0);
    ofs << "P6\n" << plot_size << " " << plot_size << "\n255\n";
    for (unsigned int y = 0; y < plot_size; y++) {
        for (unsigned int x = 0; x < plot_size; x++) {
            unsigned char c = static_cast<unsigned char>(200);
            Vec2d coordinate = Vec2d(
                ((static_cast<double>(x) / static_cast<double>(plot_size) - 0.5) * 2.0) * starshaped.max_radial_span ,
                ((static_cast<double>(y) / static_cast<double>(plot_size) - 0.5) * 2.0) * starshaped.max_radial_span);
            if (coordinate.x == 0.0 && coordinate.y == 0.0) {
                c = static_cast<unsigned char>(0);
            } else {
                if (get_distance(coordinate, origin, starshaped) < 1.0) {
                   c = static_cast<unsigned char>(0);
                }
                if (coordinate.length() > starshaped.max_radial_span) {
                    c = static_cast<unsigned char>(255);
                }
            }
            ofs << c << c << c;
        }
    }
    ofs.close();
}

void starshaped_contour(cairo_t *ctx, unsigned int plot_size, starshaped_struct& starshaped) {
    // plot starshaped set boundary (approximated by a polygonal line)
    double angle_step = 0.02;
    bool first = false;
    for (double angle = 0.0; angle < 2 * M_PI; angle+= angle_step) {
        Vec2d direction = Vec2d(cos(angle), sin(angle));
        double radial_span = get_interpolated_radial_span(direction, starshaped);
        radial_span = radial_span / starshaped.max_radial_span * starshaped_distance_cairo_length_factor;  // normalize to [0,starshaped_distance_cairo_length_factor]
        direction.x = direction.x * radial_span;
        direction.y = direction.y * radial_span;
        if (first) {
            cairo_move_to(ctx, (direction.x + 1) * plot_size * 0.5, (direction.y + 1) * plot_size * 0.5);
            first = false;
        } else {
            cairo_line_to(ctx, (direction.x + 1) * plot_size * 0.5, (direction.y + 1) * plot_size * 0.5);
        }
    }
}

void plot_starshaped_distance_cairo(unsigned int plot_size, starshaped_struct& starshaped, std::string& filename, bool is_starshaped_distance, bool plot_symmetry_axes, bool plot_interpolation_points, bool plot_interpolated_starshaped_set) {
    // Creating a cairo PDF Surface
    cairo_surface_t *csurface = cairo_pdf_surface_create(filename.c_str(), plot_size, plot_size);
    // Creating a cairo context
    cairo_t *ctx = cairo_create(csurface);
    // plot max radial span
    cairo_save(ctx);
    cairo_set_source_rgb(ctx, 0, 0, 0);
    cairo_set_line_width(ctx, plot_size * 0.005);
    static const double dashed1[] = {plot_size * 0.01, plot_size * 0.007};
    static int len1 = sizeof(dashed1) / sizeof(dashed1[0]);
    cairo_set_dash(ctx, dashed1, len1, 0);
    cairo_arc(ctx, plot_size * 0.5, plot_size * 0.5, plot_size * 0.5 * starshaped_distance_cairo_length_factor, 0, 2*M_PI);
    cairo_stroke(ctx);
    if (plot_symmetry_axes) {
        // show symmetry axis (if any)
        for (unsigned int i = 0; i < starshaped.symmetry_degree; i++) {
            double angle = static_cast<double>(i) / static_cast<double>(starshaped.symmetry_degree) * 2.0 * M_PI;
            cairo_move_to(ctx, plot_size * 0.5, plot_size * 0.5);
            cairo_line_to(ctx, plot_size * (0.5 + cos(angle) * 0.5 * starshaped_distance_cairo_length_factor), plot_size * (0.5 + sin(angle) * 0.5 * starshaped_distance_cairo_length_factor));
            cairo_stroke(ctx);
        }
    }
    cairo_restore(ctx);

    if (plot_interpolated_starshaped_set) {
        if (!plot_interpolation_points) {
            cairo_set_line_width(ctx, 0);
            starshaped_contour(ctx, plot_size, starshaped);
            if (is_starshaped_distance) {
                cairo_set_source_rgb(ctx, color_starshaped_distance[0], color_starshaped_distance[1], color_starshaped_distance[2]);
            } else {
                cairo_set_source_rgb(ctx, color_euclidean_neighbourhood[0], color_euclidean_neighbourhood[1], color_euclidean_neighbourhood[2]);
            }
        }
        cairo_fill(ctx);
        cairo_close_path(ctx);
        cairo_stroke(ctx);
        cairo_set_line_width(ctx, plot_size * 0.005);
        cairo_set_source_rgb(ctx, 0.0, 0.0, 0.0);
        starshaped_contour(ctx, plot_size, starshaped);
        cairo_close_path(ctx);
        cairo_stroke(ctx);
    }

    if (plot_interpolation_points) {
        // plot interpolating values
        for (unsigned int i = 0; i < starshaped.radial_spans.size(); i++) {
            double angle = static_cast<double>(i) / static_cast<double>(starshaped.radial_spans.size()) * (2.0 * M_PI);
            double radial_span = starshaped.radial_spans[i];
            Vec2d direction = Vec2d(cos(angle), sin(angle));
            radial_span = radial_span / starshaped.max_radial_span * starshaped_distance_cairo_length_factor;  // normalize to [0,starshaped_distance_cairo_length_factor]
            direction.x = direction.x * radial_span;
            direction.y = direction.y * radial_span;
            if (!plot_interpolated_starshaped_set) {
                cairo_set_line_width(ctx, plot_size * 0.001);
                cairo_set_source_rgb(ctx, 0.0, 0.0, 0.0);
                cairo_move_to(ctx, plot_size * 0.5, plot_size * 0.5);
                cairo_line_to(ctx, (direction.x + 1) * plot_size * 0.5, (direction.y + 1) * plot_size * 0.5);
                cairo_stroke(ctx);
            }
            cairo_set_source_rgb(ctx, 1.0, 1.0, 1.0);
            cairo_arc(ctx, (direction.x + 1) * plot_size * 0.5, (direction.y + 1) * plot_size * 0.5, plot_size * 0.007, 0, 2*M_PI);
            cairo_fill(ctx);
            cairo_set_line_width(ctx, plot_size * 0.004);
            cairo_set_source_rgb(ctx, 0.0, 0.0, 0.0);
            cairo_arc(ctx, (direction.x + 1) * plot_size * 0.5, (direction.y + 1) * plot_size * 0.5, plot_size * 0.007, 0, 2*M_PI);
            cairo_stroke(ctx);
        }
    }

    // plot origin point
    cairo_set_source_rgb(ctx, 0, 0, 0);
    cairo_set_line_width(ctx, 0.0);
    cairo_arc(ctx, plot_size * 0.5, plot_size * 0.5, plot_size * 0.007, 0, 2*M_PI);
    cairo_fill(ctx);
    cairo_stroke(ctx);

    // Destroying cairo context
    cairo_destroy(ctx);
    cairo_surface_flush(csurface);
    // Destroying PDF surface
    cairo_surface_destroy(csurface);
    cairo_debug_reset_static_data();
}

void plot_sites(cairo_t *ctx, std::vector<Vec2d>& sites, unsigned int dim_x, unsigned int dim_y) {
    cairo_set_source_rgb(ctx, color_point_process[0], color_point_process[1], color_point_process[2]);
    cairo_set_line_width(ctx, 0.0);
    for (unsigned int i = 0; i < sites.size(); i++) {
        cairo_arc(ctx, dim_x * sites[i].x, dim_y * sites[i].y, dim_x * 0.007, 0, 2*M_PI);
        cairo_fill(ctx);
    }
}

void plot_sites_pdf(std::vector<Vec2d>& sites, unsigned int dim_x, unsigned int dim_y, std::string& filename) {
    cairo_surface_t *csurface = cairo_pdf_surface_create(filename.c_str(), dim_x, dim_y);
    cairo_t *ctx = cairo_create(csurface);
    plot_sites(ctx, sites, dim_x, dim_y);
    cairo_destroy(ctx);
    cairo_surface_flush(csurface);
    cairo_surface_destroy(csurface);
    cairo_debug_reset_static_data();
}

void plot_porous_material_and_sites_png(std::vector<Vec2d>& sites, unsigned int dim_x, unsigned int dim_y, std::string& filename, const std::vector<unsigned int>& image, bool is_porous_material_growing, bool is_plot_sites) {
    cairo_surface_t *csurface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, dim_x, dim_y);
    cairo_t *ctx = cairo_create(csurface);
    // draw background
    if (is_porous_material_growing) {
        cairo_set_source_rgb(ctx, color_porous_material_growing[0], color_porous_material_growing[1], color_porous_material_growing[2]);
    } else {
        cairo_set_source_rgb(ctx, color_porous_material[0], color_porous_material[1], color_porous_material[2]);
    }
    cairo_rectangle(ctx, 0, 0, dim_x, dim_y);
    cairo_fill(ctx);
    // draw discrete cells
     if (is_porous_material_growing) {
        cairo_set_source_rgb(ctx, color_cell_growing[0], color_cell_growing[1], color_cell_growing[2]);
    } else {
        cairo_set_source_rgb(ctx, 1.0, 1.0, 1.0);
    }
    for (unsigned int y = 0; y < dim_y; y++) {
        for (unsigned int x = 0; x < dim_x ; x++) {
            bool fill = false;
            if (is_porous_material_growing) {
                fill = (image[x + y * dim_x] != std::numeric_limits<unsigned int>::max());
            } else {
                fill = (image[x + y * dim_x] == void_phase);
            }
            if (fill) {
                cairo_rectangle(ctx, x, y, 1, 1);
                cairo_fill(ctx);
            }
        }
    }
    if (is_plot_sites) {
        plot_sites(ctx, sites, dim_x, dim_y);
    }
    cairo_destroy(ctx);
    cairo_surface_flush(csurface);
    cairo_surface_write_to_png(csurface, filename.c_str());
    cairo_surface_destroy(csurface);
    cairo_debug_reset_static_data();
}

void save_porous_material_ppm(std::string& filename, const std::vector<unsigned int>& image_porous_material, unsigned int dim_x, unsigned int dim_y) {
    std::ofstream ofs;
    ofs.open(filename.c_str(), std::ios::binary);
    ofs << "P6\n" << dim_x << " " << dim_y << "\n255\n";
    for (unsigned int y = 0; y < dim_y; y++) {
        for (unsigned int x = 0; x < dim_x ; x++) {
            unsigned char c = static_cast<unsigned char>(0);
            if (image_porous_material[x + y * dim_x] == solid_phase) {
                c = static_cast<unsigned char>(255);
            }
            ofs << c << c << c;
        }
    }
    ofs.close();
}

#endif

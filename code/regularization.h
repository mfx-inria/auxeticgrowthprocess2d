#ifndef REGULARIZATION_H
#define REGULARIZATION_H

#include <queue>
#include <vector>

// convert from (c, oc) pair, e.g. x and osx, to absolute position
int abs_coord(int c, int oc, int s) {
    return c - oc * s;
}

void get_abs_coords(int &abs_coord_x, int& abs_coord_y, const GrowthInfo& growth_info, unsigned int dim_x, unsigned int dim_y) {
    abs_coord_x = abs_coord(growth_info.x, growth_info.osx, static_cast<int>(dim_x));
    abs_coord_y = abs_coord(growth_info.y, growth_info.osy, static_cast<int>(dim_y));
}

void regularize_cells(std::vector<std::vector<GrowthInfo> >& Cw, std::vector<unsigned int>& image_porous_material, unsigned int dim_x, unsigned int dim_y, unsigned int& num_void_pixels_regularization) {
    num_void_pixels_regularization = 0;
    std::fill(image_porous_material.begin(), image_porous_material.end(), solid_phase);
    #pragma omp parallel for
    for (unsigned int i = 0; i < Cw.size(); i++) {
        std::vector<GrowthInfo> Cwi = Cw[i];
        // retrieve the minimum and maximum coordinate of the i-th cell
        GrowthInfo max_growth_info, min_growth_info;
        max_growth_info.x = max_growth_info.y = -1;
        max_growth_info.osx = max_growth_info.osy = min_growth_info.osx = min_growth_info.osy = 0;
        min_growth_info.x =  dim_x + 1;
        min_growth_info.y = dim_y + 1;
        for (unsigned int c = 0; c < Cwi.size(); c++) {
            GrowthInfo Cwi_c = Cwi[c];
            int abs_coord_max_x, abs_coord_min_x, abs_coord_max_y, abs_coord_min_y, abs_coord_x, abs_coord_y;
            get_abs_coords(abs_coord_max_x, abs_coord_max_y, max_growth_info, dim_x, dim_y);
            get_abs_coords(abs_coord_min_x, abs_coord_min_y, min_growth_info, dim_x, dim_y);
            get_abs_coords(abs_coord_x, abs_coord_y, Cwi_c, dim_x, dim_y);
            if (abs_coord_x > abs_coord_max_x) {max_growth_info.x = Cwi_c.x; max_growth_info.osx = Cwi_c.osx;}
            if (abs_coord_y > abs_coord_max_y) {max_growth_info.y = Cwi_c.y; max_growth_info.osy = Cwi_c.osy;}
            if (abs_coord_x < abs_coord_min_x) {min_growth_info.x = Cwi_c.x;  min_growth_info.osx = Cwi_c.osx;}
            if (abs_coord_y < abs_coord_min_y) { min_growth_info.y = Cwi_c.y; min_growth_info.osy = Cwi_c.osy;}
        }

        // flood filling for hole filling, with an outer border of one pixel size
        int abs_coord_max_x, abs_coord_min_x, abs_coord_max_y, abs_coord_min_y;
        get_abs_coords(abs_coord_max_x, abs_coord_max_y, max_growth_info, dim_x, dim_y);
        get_abs_coords(abs_coord_min_x, abs_coord_min_y, min_growth_info, dim_x, dim_y);

        int cell_size_x = abs_coord_max_x + 1 + 2 - abs_coord_min_x;
        int cell_size_y = abs_coord_max_y + 1 + 2 - abs_coord_min_y;
        if (cell_size_x > 0 && cell_size_y > 0) {
            std::vector<unsigned int> image_cell(cell_size_x * cell_size_y);
            std::fill(image_cell.begin(), image_cell.end(), 1);
            for (unsigned int c = 0; c < Cwi.size(); c++) {
                GrowthInfo Cwi_c = Cwi[c];
                int cx = abs_coord(Cwi_c.x, Cwi_c.osx, static_cast<int>(dim_x)) - abs_coord_min_x + 1;
                int cy = abs_coord(Cwi_c.y, Cwi_c.osy, static_cast<int>(dim_y)) - abs_coord_min_y + 1;
                image_cell[cx + cy * cell_size_x] = 0;
            }
            std::queue<unsigned int> to_visit;
            to_visit.push(0);
            std::vector<unsigned int> image_visited(cell_size_x * cell_size_y);
            std::fill(image_visited.begin(), image_visited.end(), 0);  // flood fill from the boundary, point (0,0)
            while (!to_visit.empty()) {
                unsigned int pv = to_visit.front();
                unsigned int xv = pv % cell_size_x;
                unsigned int yv = pv / cell_size_x;
                to_visit.pop();
                if (image_visited[pv] == 0) {
                    image_visited[pv] = 1;
                    for (unsigned int k = 0; k < 4; k++) {
                        int xn = xv + neigh_indices[k * 2];
                        int yn = yv + neigh_indices[k * 2 + 1];
                        if ((xn >= 0 && xn < cell_size_x) && (yn >= 0 && yn < cell_size_y)) {
                            unsigned int pn = xn + yn * cell_size_x;
                            if (image_cell[pn] == 1) {
                                to_visit.push(pn);
                            }
                        }
                    }
                }
            }
            // copy cell back in the image (take into account border)
            for (int x = 1; x < cell_size_x - 1; x++) {
                for (int y = 1; y < cell_size_y - 1; y++) {
                    if (image_visited[x + cell_size_x * y] == 0) {
                        int xc = x + abs_coord_min_x - 1;
                        int osx = 0, osy = 0;
                        funcmod(xc, osx, dim_x);
                        int yc = y + abs_coord_min_y - 1;
                        funcmod(yc, osy, dim_y);
                        image_porous_material[xc + dim_x * yc] = void_phase;
                        #pragma omp atomic
                        num_void_pixels_regularization++;
                    }
                }
            }
        }
    }
}

#endif

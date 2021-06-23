#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>


struct parameters_struct {
    std::string name;
    unsigned int image_size;
    starshaped_struct starshaped;
    double max_growth_length;
    double max_radial_span, max_length_growth_radial_span;  // input max radial spans (not necessarily being the max among the input radial spans)
    starshaped_struct max_growth_length_starshaped;
    PointProcessType point_process;
    SymmetryType symmetry_type;
    unsigned int symmetry_degree, random_seed, num_growth_plots;
    double rsa_max_dist_rel_pixel, num_sites_per_pixel;
    std::string filename_points;
    bool plot_symmetry_axes, save_non_regularized, plot_starshaped_ppm, plot_starshaped_pdf, plot_sites_pdf, save_porous_material_ppm, plot_porous_material_sites_png, save_results_txt, plot_interpolation_points, plot_interpolated_starshaped_set, plot_sites_png;
};

void string_list_to_vector(std::vector<std::string>& l, std::string s, std::string delimiter) {
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        l.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    l.push_back(s);
}

void read_radial_spans(const std::string& value, std::vector<double>& radial_spans) {
    std::vector<std::string> radial_spans_str;
    string_list_to_vector(radial_spans_str, value, ",");
    for (unsigned int i = 0; i < radial_spans_str.size(); i++) {
        double radial_span = static_cast<double>(atof(radial_spans_str[i].c_str()));
        radial_spans.push_back(radial_span);
    }
}

bool to_bool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b = false;
    is >> std::boolalpha >> b;
    return b;
}


bool initialize_parameters(const std::string& filename, parameters_struct& parameters) {
    parameters.image_size = 0;
    parameters.max_growth_length = 2.0;
    parameters.symmetry_degree = 0;
    parameters.max_radial_span = parameters.max_length_growth_radial_span = 0.0;
    parameters.random_seed = 2;
    parameters.num_sites_per_pixel = 0.0;
    parameters.num_growth_plots = 0;
    parameters.rsa_max_dist_rel_pixel = 0.0;
    parameters.save_non_regularized = false;
    parameters.plot_interpolation_points = false;
    parameters.plot_interpolated_starshaped_set = true;
    parameters.plot_sites_png = true;
    parameters.plot_starshaped_ppm = parameters.save_porous_material_ppm = parameters.save_results_txt = true;
    parameters.plot_starshaped_pdf = parameters.plot_porous_material_sites_png = parameters.plot_sites_pdf = parameters.plot_symmetry_axes = false;
    std::ifstream parameters_file(filename.c_str());
    std::string line;
    std::cout << "* Parameters" << std::endl;
    while (std::getline(parameters_file, line)) {
        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, '=')) {
            std::string value;
            if (std::getline(is_line, value)) {
                if (key.compare("image_size") == 0) parameters.image_size = static_cast<unsigned int>(atoi(value.c_str()));
                else if (key.compare("max_growth_length") == 0) parameters.max_growth_length = static_cast<double>(atof(value.c_str()));
                else if (key.compare("max_radial_span") == 0) parameters.max_radial_span = static_cast<double>(atof(value.c_str()));
                else if (key.compare("max_length_growth_radial_span") == 0) parameters.max_length_growth_radial_span = static_cast<double>(atof(value.c_str()));
                else if (key.compare("symmetry_degree") == 0) parameters.symmetry_degree = static_cast<unsigned int>(atoi(value.c_str()));
                else if (key.compare("random_seed") == 0) parameters.random_seed = static_cast<unsigned int>(atoi(value.c_str()));
                else if (key.compare("num_sites_per_pixel") == 0) parameters.num_sites_per_pixel = static_cast<double>(atof(value.c_str()));
                else if (key.compare("num_growth_plots") == 0) parameters.num_growth_plots = static_cast<unsigned int>(atoi(value.c_str()));
                else if (key.compare("rsa_max_dist_rel_pixel") == 0) parameters.rsa_max_dist_rel_pixel = static_cast<double>(atof(value.c_str()));
                else if (key.compare("name") == 0) parameters.name = value;
                else if (key.compare("point_process") == 0) {
                    if (value.compare("RSA") == 0) parameters.point_process = RSA;
                    else if (value.compare("Random") == 0) parameters.point_process = Random;
                    else if (value.compare("File") == 0) parameters.point_process = File;
                } else if (key.compare("radial_spans") == 0) {
                    read_radial_spans(value, parameters.starshaped.radial_spans);
                } else if (key.compare("max_growth_length_radial_spans") == 0) {
                    read_radial_spans(value, parameters.max_growth_length_starshaped.radial_spans);
                } else if (key.compare("interpolation_type") == 0) {
                    if (value.compare("PolarPiecewise") == 0) {
                        parameters.starshaped.interpolation_type = parameters.max_growth_length_starshaped.interpolation_type = PolarPiecewise;
                    } else if (value.compare("PolarLinear") == 0) {
                        parameters.starshaped.interpolation_type = parameters.max_growth_length_starshaped.interpolation_type =  PolarLinear;
                    } else if (value.compare("PolarCubic") == 0) {
                        parameters.starshaped.interpolation_type = parameters.max_growth_length_starshaped.interpolation_type = PolarCubic;
                    } else if (value.compare("Polygonal") == 0) {
                        parameters.starshaped.interpolation_type = parameters.max_growth_length_starshaped.interpolation_type = Polygonal;
                    }
                } else if (key.compare("symmetry_type") == 0) {
                    if (value.compare("NoSymmetry") == 0) parameters.symmetry_type = NoSymmetry;
                    else if (value.compare("RotationalSymmetry") == 0) parameters.symmetry_type = RotationalSymmetry;
                    else if (value.compare("ReflectionalSymmetry") == 0) parameters.symmetry_type = ReflectionalSymmetry;
                } else if (key.compare("plot_starshaped_ppm") == 0) {parameters.plot_starshaped_ppm = to_bool(value);}
                else if (key.compare("plot_starshaped_pdf") == 0) {parameters.plot_starshaped_pdf = to_bool(value);}
                else if (key.compare("plot_sites_pdf") == 0) {parameters.plot_sites_pdf = to_bool(value);}
                else if (key.compare("save_non_regularized") == 0) {parameters.save_non_regularized = to_bool(value);}
                else if (key.compare("plot_symmetry_axes") == 0) {parameters.plot_symmetry_axes = to_bool(value);}
                else if (key.compare("save_porous_material_ppm") == 0) {parameters.save_porous_material_ppm = to_bool(value);}
                else if (key.compare("save_results_txt") == 0) {parameters.save_results_txt = to_bool(value);}
                else if (key.compare("plot_sites_png") == 0) {parameters.plot_sites_png = to_bool(value);}
                else if (key.compare("plot_interpolation_points") == 0) {parameters.plot_interpolation_points = to_bool(value);}
                else if (key.compare("plot_interpolated_starshaped_set") == 0) {parameters.plot_interpolated_starshaped_set = to_bool(value);}
                else if (key.compare("plot_porous_material_sites_png") == 0) {parameters.plot_porous_material_sites_png = to_bool(value);}
                else if (key.compare("filename_points") == 0) {parameters.filename_points = value;}
                std::cout << "\t" << line << std::endl;
            }
        }
    }
    parameters_file.close();
    bool valid_parameters = true;
    const double min_growth_length = sqrt(2.0);
    if (parameters.image_size == 0) {
        std::cerr << "[error] image_size must be greater than zero" << std::endl;
        valid_parameters = false;
    }
    if (parameters.name.empty()) {
        std::cerr << "[error] name not provided" << std::endl;
        valid_parameters = false;
    }
    if (parameters.starshaped.radial_spans.empty()) {
        std::cerr << "[error] empty radial spans" << std::endl;
        valid_parameters = false;
    }
    for (unsigned int i = 0; i < parameters.starshaped.radial_spans.size(); i++) {
        if (parameters.starshaped.radial_spans[i] <= 0.0) {
            std::cerr << "[error] nonpositive radial span " << parameters.starshaped.radial_spans[i] << std::endl;
            valid_parameters = false;
        }
    }
    if (!parameters.max_growth_length_starshaped.radial_spans.empty()) {
        if (parameters.max_growth_length_starshaped.radial_spans.size() != parameters.starshaped.radial_spans.size()) {
            std::cerr << "[error] not equal number of radial spans (max growth length)" << std::endl;
            valid_parameters = false;
        }
        for (unsigned int i = 0; i < parameters.max_growth_length_starshaped.radial_spans.size(); i++) {
            if (parameters.max_growth_length_starshaped.radial_spans[i] < min_growth_length) {
                std::cerr << "[error] max_growth_length (radial spans) < " <<  min_growth_length << std::endl;
                valid_parameters = false;
            }
        }
    }
    if (parameters.max_growth_length < min_growth_length) {
        std::cerr << "[error] max_growth_length < " << min_growth_length << std::endl;
        valid_parameters = false;
    }
    if (parameters.random_seed < 2) {
        std::cerr << "[error] random_seed < 2" << std::endl;
        valid_parameters = false;
    }
    if (parameters.point_process == Random) {
        if (parameters.num_sites_per_pixel <= 0.0) {
            std::cerr << "[error] nonpositive num_sites_per_pixel" << std::endl;
            valid_parameters = false;
        }
    } else if (parameters.point_process == RSA) {
        if (parameters.rsa_max_dist_rel_pixel <= 0.0) {
            std::cerr << "[error] nonpositive rsa_max_dist_rel_pixel" << std::endl;
            valid_parameters = false;
        }
    }
    if (parameters.symmetry_type != NoSymmetry) {
        if (parameters.symmetry_degree < 1) {
            std::cerr << "[error] symmetry_degree < 1" << std::endl;
            valid_parameters = false;
        }
        if (parameters.symmetry_type == ReflectionalSymmetry) {
            if ((parameters.symmetry_degree % 2) != 0) {
                std::cerr << "[error] (symmetry_degree % 2) != 0" << std::endl;
                valid_parameters = false;
            }
        }
    }
    return valid_parameters;
}

#endif

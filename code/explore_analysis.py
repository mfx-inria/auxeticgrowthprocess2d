import sys
import os
import numpy

import utils


def explore_analysis_main(filename):
    # read data
    parameters = utils.read_parameters_explore(filename)
    folder_name = parameters["name"]
    all_radial_variables = utils.get_radial_variables_random(parameters)
    deviation_isotropy_set_dict_image_size = {}
    porosity_coefficient_variation_set_dict_image_size = {}
    young_set_dict_image_size = {}
    poisson_set_dict_image_size = {}
    range_image_size = utils.get_discrete_set(
        int(parameters["image_size_ini"]), int(parameters["image_size_end"]),
        int(parameters["num_samples_image_size"]))
    range_image_size = map(int, range_image_size)  # list to int
    porosity_regularized_set_all = []
    porosity_set_all = []
    for image_size in range_image_size:
        print("*** image_size=" + str(image_size))
        num_computed = 0
        deviation_isotropy_set = []
        porosity_coefficient_variation_set = []
        young_set = []
        poisson_set = []
        parameters["image_size"] = str(image_size)
        id_variables = 0
        for variables in all_radial_variables:
            # print("* Test id = " + str(id_variables))
            name = folder_name + "image_size_" + str(image_size) + "_" + str(id_variables)
            parameters["name"] = name
            utils.parse_porous_material_from_variables(parameters, variables)
            (elasticity_tensor_set, porosity_set, porosity_regularized_set, num_connected_components_set, num_connected_components_regularized_set, num_sites_set) = utils.parse_compute_porous_material_stochastic_homogenize_2d(parameters)
            for c in range(len(num_connected_components_regularized_set)):
                if num_connected_components_regularized_set[c] != 1:
                    print("Error: num_connected_components_regularized_set!= 1 for test id = " + str(str(id_variables)) + " and sample " + str(c))
            if (len(elasticity_tensor_set) > 0):
                isotropic_elasticity_tensor = utils.get_closest_isotropic_elasticity_tensor_frobenius(numpy.mean(elasticity_tensor_set, axis=0))
                poisson = utils.get_poissons_ratio(isotropic_elasticity_tensor)
                young = utils.get_young_modulus(isotropic_elasticity_tensor)
                poisson_set.append(poisson)
                young_set.append(young)
            porosity_set_all.extend(porosity_set)
            porosity_regularized_set_all.extend(porosity_regularized_set)
            if len(elasticity_tensor_set) == int(parameters["num_samples_montecarlo"]):
                deviation_isotropy_set.append(utils.get_deviation_isotropy(numpy.mean(elasticity_tensor_set, axis=0)))
                num_computed += 1
                porosity_coefficient_variation = utils.coefficient_variation(porosity_regularized_set)
                porosity_coefficient_variation_set.append(porosity_coefficient_variation)
            id_variables += 1
        print("\t- Num computed = " + str(num_computed))
        if len(deviation_isotropy_set) > 0:
            print("\t- Maximum deviation isotropy of " + str(numpy.max(deviation_isotropy_set)) + " with id " + str(numpy.argmax(deviation_isotropy_set)))
            print("\t- Minimum deviation isotropy of " + str(numpy.min(deviation_isotropy_set)) + " with id " + str(numpy.argmin(deviation_isotropy_set)))
            print("\t- Mean deviation isotropy of " + str(numpy.mean(deviation_isotropy_set)))
            deviation_isotropy_set_dict_image_size[int(image_size)] = deviation_isotropy_set
        if len(porosity_coefficient_variation_set) > 0:
            print("\t- Maximum porosity coefficient variation of " + str(numpy.max(porosity_coefficient_variation_set)) + " with id " + str(numpy.argmax(porosity_coefficient_variation_set)))
            print("\t- Minimum porosity coefficient variation of " + str(numpy.min(porosity_coefficient_variation_set)) + " with id " + str(numpy.argmin(porosity_coefficient_variation_set)))
            print("\t- Mean porosity coefficient variation of " + str(numpy.mean(porosity_coefficient_variation_set)))
            porosity_coefficient_variation_set_dict_image_size[int(image_size)] = porosity_coefficient_variation_set
        if len(poisson_set) > 0:
            print("\t- Maximum Poisson ratio of " + str(numpy.max(poisson_set)) + " with id " + str(numpy.argmax(poisson_set)))
            print("\t- Minimum Poisson ratio of " + str(numpy.min(poisson_set)) + " with id " + str(numpy.argmin(poisson_set)))
            poisson_set_dict_image_size[int(image_size)] = poisson_set
        if len(poisson_set) > 0:
            print("\t- Maximum Young modulus of " + str(numpy.max(young_set)) + " with id " + str(numpy.argmax(young_set)))
            print("\t- Minimum Young modulus of " + str(numpy.min(young_set)) + " with id " + str(numpy.argmin(young_set)))
            young_set_dict_image_size[int(image_size)] = young_set

    if len(porosity_regularized_set_all) > 0:
        print("- Range of porosity (regularized) = " + str(numpy.min(porosity_regularized_set_all)) + "," + str(numpy.max(porosity_regularized_set_all)))
        porosity_set_all = numpy.array(porosity_set_all)
        porosity_regularized_set_all = numpy.array(porosity_regularized_set_all)
        dev_porosity_regularized = numpy.abs(numpy.divide(porosity_set_all - porosity_regularized_set_all, porosity_set_all))
        max_dev_porosity_regularized = numpy.max(dev_porosity_regularized)
        print("- Maximum deviation porosity regularized = " + str(max_dev_porosity_regularized))
        mean_dev_porosity_regularized = numpy.mean(dev_porosity_regularized)
        print("- Mean deviation porosity regularized = " + str(mean_dev_porosity_regularized))

    parameters_single_optimize = utils.read_parameters("parameters/global_single_optimize.txt")
    target_image_size = int(parameters_single_optimize["image_size"])
    print("* target_image_size = " + str(target_image_size))
    assert target_image_size in porosity_coefficient_variation_set_dict_image_size
    assert target_image_size in deviation_isotropy_set_dict_image_size

    target_max_coeff_var_porosity = numpy.max(porosity_coefficient_variation_set_dict_image_size[target_image_size])
    print("* Target maxium coefficient variation porosity = " + str(target_max_coeff_var_porosity))
    target_max_dev_iso = numpy.max(deviation_isotropy_set_dict_image_size[target_image_size])
    print("* Target maxium deviation isotropy = " + str(target_max_dev_iso))

    with open(os.path.join(folder_name, "latex.txt"), 'w') as f:
        utils.write_latex_command(f, parameters["point_process"] + "maxdevpor", max_dev_porosity_regularized)
        utils.write_latex_command(f, parameters["point_process"] + "meandevpor", mean_dev_porosity_regularized)
        utils.write_latex_command(f, parameters["point_process"] + "targetmaxcoeffvarpor", target_max_coeff_var_porosity)
        utils.write_latex_command(f, parameters["point_process"] + "targetmaxdeviso", target_max_dev_iso)

    return (folder_name, deviation_isotropy_set_dict_image_size, porosity_coefficient_variation_set_dict_image_size)


if __name__ == "__main__":
    assert(len(sys.argv) == 2)
    explore_analysis_main(sys.argv[1])

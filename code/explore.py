import sys

import utils

if __name__ == "__main__":
    parameters = utils.read_parameters_explore(sys.argv[1])
    folder_name = parameters["name"]
    all_radial_variables = utils.get_radial_variables_random(parameters)
    range_image_size = utils.get_discrete_set(int(parameters["image_size_ini"]), int(parameters["image_size_end"]), int(parameters["num_samples_image_size"]))
    range_image_size = map(int, range_image_size)  # list to int
    for image_size in range_image_size:
        parameters["image_size"] = str(image_size)
        id_variables = 0
        for variables in all_radial_variables:
            print("* Test id = " + str(id_variables))
            name = folder_name + "image_size_" + str(image_size) + "_" + str(id_variables)
            parameters["name"] = name
            utils.parse_porous_material_from_variables(parameters, variables)
            utils.compute_porous_material_stochastic_homogenize_2d(parameters)
            id_variables += 1

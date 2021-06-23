import sys
import os
import numpy

import utils


def prime_factors(n):
    f = 2
    factors = []
    while f * f <= n:
        if n % f:
            f += 1
        else:
            n //= f
            factors.append(f)
    if n > 1:
        factors.append(n)
    return factors


def get_image_sizes(pixel_size_ini, pixel_size_end, num_samples_timing, num_factors_two=4):
    # we consider images with size with specific number of factors of two in order to lessen the impact
    # of data structure alignment
    image_sizes = []
    for i in range(num_samples_timing):
        image_size_target = int(numpy.sqrt(pixel_size_ini + float(i) / float(num_samples_timing - 1) * (pixel_size_end - pixel_size_ini)))
        while (prime_factors(image_size_target).count(2) != num_factors_two):
            image_size_target = image_size_target + 1
        image_sizes.append(image_size_target)
    return image_sizes


def parse_parameters_timings(filename_parameters):
    parameters = utils.read_parameters_timings(filename_parameters)
    image_size_ini = int(parameters["image_size_ini"])
    image_size_end = int(parameters["image_size_end"])
    num_samples_timing = int(parameters["num_samples_timing"])
    folder_name = os.path.join(os.path.dirname(parameters["name"]), "timings")
    pixel_size_ini = image_size_ini ** 2
    pixel_size_end = image_size_end ** 2
    parameters["plot_starshaped_pdf"] = "false"
    parameters["plot_porous_material_sites_png"] = "false"
    parameters["save_porous_material_ppm"] = "false"
    range_growth_length = None
    if "min_length_growth_radial_span" in parameters and "max_length_growth_radial_span" in parameters:
        range_growth_length = [parameters["min_length_growth_radial_span"], parameters["max_length_growth_radial_span"]]
    elif "max_growth_length" in parameters:
        range_growth_length = [parameters["max_growth_length"]]
    return (num_samples_timing, pixel_size_ini, pixel_size_end, folder_name, range_growth_length, parameters)


if __name__ == "__main__":
    assert len(sys.argv) == 2
    (num_samples_timing, pixel_size_ini, pixel_size_end, folder_name, range_growth_length, parameters) = parse_parameters_timings(sys.argv[1])
    image_sizes = get_image_sizes(pixel_size_ini, pixel_size_end, num_samples_timing)
    for i, image_size in enumerate(image_sizes):
        print(image_size)
        parameters["image_size"] = str(image_size)
        parameters["name"] = os.path.join(folder_name, str(i))
        print("* Test id = " + str(i) + " image_size=" + parameters["image_size"])
        utils.compute_porous_material(parameters)

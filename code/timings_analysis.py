import numpy
import os

import utils
import timings


def format_number(num):
    return numpy.format_float_positional(num, precision=3, fractional=False, trim='-')


def timings_analysis_main(filename, file_latex):
    print("*** Timings analysis " + filename)
    (num_samples_timing, pixel_size_ini, pixel_size_end, folder_name, range_growth_length, parameters) = timings.parse_parameters_timings(filename)
    pixels_million_set = []
    time_growth_process_min_set = []
    time_cell_regularization_min_set = []
    num_iterations_growth_set = []
    avg_size_queue_set = []
    avg_iterations_loop_discard_set = []
    image_sizes = timings.get_image_sizes(pixel_size_ini, pixel_size_end, num_samples_timing)
    for i, image_size in enumerate(image_sizes):
        parameters["image_size"] = str(image_size)
        parameters["name"] = os.path.join(folder_name, str(i))
        filename_results = parameters["name"] + ".res"
        if os.path.isfile(filename_results):
            results = utils.read_parameters(filename_results)
            time_growth_process_ms = int(results["time_growth_process_ms"])
            time_cell_regularization_ms = int(results["time_cell_regularization_ms"])
            pixels_million_set.append(float(image_size * image_size) / 1000000.0)
            time_growth_process_min_set.append(float(time_growth_process_ms) / (1000.0 * 60.0))
            time_cell_regularization_min_set.append(float(time_cell_regularization_ms) / (1000.0 * 60.0))
            num_iterations_growth_set.append(float(results["num_iterations_growth"]) * 100 / float(image_size * image_size))
            avg_size_queue_set.append(float(results["avg_size_queue"]) * 100.0 / float(image_size * image_size))
            avg_iterations_loop_discard_set.append(float(results["avg_iterations_loop_discard"]))
    dirname = os.path.dirname(filename)
    name_latex = (dirname.split("/")[-2] + dirname.split("/")[-1]).replace("_", "")
    num_iterations_growth_mean = numpy.mean(num_iterations_growth_set)
    avg_size_queue_mean = numpy.mean(avg_size_queue_set)
    avg_iterations_loop_discard_mean = numpy.mean(avg_iterations_loop_discard_set)
    print("* Mean num_iterations_growth (%) = " + str(num_iterations_growth_mean))
    print("\t-Coefficient variation num_iterations_growth = " + str(utils.coefficient_variation(num_iterations_growth_set)))
    print("* Mean avg_size_queue (%) = " + str(avg_size_queue_mean))
    print("\t-Coefficient variation avg_size_queue = " + str(utils.coefficient_variation(avg_size_queue_set)))
    print("* Mean avg_iterations_loop_discard = " + str(avg_iterations_loop_discard_mean))
    print("\t-Coefficient variation avg_iterations_loop_discard = " + str(utils.coefficient_variation(avg_iterations_loop_discard_set)))
    pixels_million_set = numpy.array(pixels_million_set)
    time_growth_process_min_set = numpy.array(time_growth_process_min_set)
    time_cell_regularization_min_set = numpy.array(time_cell_regularization_min_set)
    time_total_min_set = time_growth_process_min_set + time_cell_regularization_min_set
    time_regularization_relative_set = numpy.divide(time_cell_regularization_min_set, time_total_min_set) * 100
    time_regularization_relative_mean = numpy.mean(time_regularization_relative_set)
    print("* Mean time_regularization_relative (%) = " + str(time_regularization_relative_mean))
    print("\t-Coefficient variation time_regularization_relative = " + str(utils.coefficient_variation(time_regularization_relative_set)))
    value_latex = format_number(time_regularization_relative_mean) + "\\%&" + format_number(num_iterations_growth_mean) + "\\%&" + format_number(avg_size_queue_mean) + "\\%&" + format_number(avg_iterations_loop_discard_mean)
    utils.write_latex_command(file_latex, name_latex, value_latex, si_format=False)
    return (folder_name, pixels_million_set, time_growth_process_min_set, range_growth_length)

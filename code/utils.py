import os
import multiprocessing
import shutil
import numpy


def write_parameters(local_parameters, filename, sort=False):
    with open(filename, "w") as f:
        if sort:
            for param in sorted(local_parameters):
                f.write(param + "=" + str(local_parameters[param]) + "\n")
        else:
            for param in local_parameters:
                f.write(param + "=" + str(local_parameters[param]) + "\n")


def read_parameters(filename, parameters=None):
    with open(filename, "r") as f:
        fields = f.read().strip().split("\n")
    if parameters is None:
        parameters = {}
    for field in fields:
        field_name = field.split("=")[0]
        field_value = field.split("=")[1]
        if field_name in parameters:
            print(str(field_name) + " in parameters")
        # assert field_name not in parameters
        parameters[field_name] = field_value
    return parameters


def write_latex_command(f, name, value, si_format=True):
    if si_format:
        text = "\\newcommand{\\" + str(name) + "}{\\SI{" + str(value) + "}{}}"
    else:
        text = "\\newcommand{\\" + str(name) + "}{" + str(value) + "}"
    f.write(text + "\n")


def read_parameters_single(filename):
    parameters = read_parameters(filename)
    parameters = read_parameters("parameters/global.txt", parameters)
    parameters = read_parameters("parameters/global_single.txt", parameters)
    parameters = read_parameters("parameters/global_single_optimize.txt", parameters)
    return parameters


def read_parameters_explore(filename):
    parameters = read_parameters(filename)
    parameters = read_parameters("parameters/global.txt", parameters)
    parameters = read_parameters("parameters/global_explore_optimize.txt", parameters)
    parameters = read_parameters("parameters/global_explore.txt", parameters)
    return parameters


def read_parameters_timings(filename):
    parameters = read_parameters(filename)
    parameters = read_parameters("parameters/global_timings.txt", parameters)
    return parameters


def read_parameters_optimize(filename):
    parameters = read_parameters(filename)
    parameters = read_parameters("parameters/global.txt", parameters)
    parameters = read_parameters("parameters/global_explore_optimize.txt", parameters)
    parameters = read_parameters("parameters/global_single_optimize.txt", parameters)
    parameters = read_parameters("parameters/global_optimize.txt", parameters)
    return parameters


def has_parameterized_length_growth(parameters):
    return ("min_length_growth_radial_span" in parameters) and ("max_length_growth_radial_span" in parameters)


def get_discrete_set(ini, end, cardinality):
    assert cardinality >= 1
    assert ini <= end
    discrete_set = [ini]
    for i in range(1, cardinality - 1):
        discrete_set.append(ini + (end - ini) / (float(cardinality) - 1.0) * float(i))
    discrete_set.append(end)
    return discrete_set


def is_lock_max_radial_span(parameters):
    if "lock_max_radial_span" in parameters:
        return parameters["lock_max_radial_span"] == "true"
    return False


def parse_porous_material_from_variables(parameters, variables, name=None):
    if has_parameterized_length_growth(parameters):
        parameters["radial_spans"] = ",".join(list(map(str, variables[:int(len(variables) / 2)])))
        parameters["max_growth_length_radial_spans"] = ",".join(list(map(str, variables[int(len(variables) / 2):])))
    else:
        parameters["radial_spans"] = ",".join(list(map(str, variables)))
    if name is not None:
        parameters["name"] = name
    if is_lock_max_radial_span(parameters):
        parameters["radial_spans"] = "1," + parameters["radial_spans"]


def get_radial_variables_random(parameters):
    numpy.random.seed(int(parameters["random_seed"]))
    min_radial_span = float(parameters["min_radial_span"])
    max_radial_span = float(parameters["max_radial_span"])
    num_radial_spans = int(parameters["num_radial_spans"])
    all_radial_variables = []
    for t in range(int(parameters["num_samples_radial_spans"])):
        radial_spans = numpy.random.uniform(
            min_radial_span, max_radial_span, num_radial_spans)
        if has_parameterized_length_growth(parameters):
            min_length_growth_radial_span = float(parameters["min_length_growth_radial_span"])
            max_length_growth_radial_span = float(parameters["max_length_growth_radial_span"])
            radial_spans_length_growth = numpy.random.uniform(
                min_length_growth_radial_span, max_length_growth_radial_span, num_radial_spans)
            radial_spans = numpy.concatenate((radial_spans, radial_spans_length_growth))
        all_radial_variables.append(radial_spans)
    return all_radial_variables


def compute_porous_material(parameters):
    filename_parameters = parameters["name"] + ".txt"
    write_parameters(parameters, filename_parameters)
    exe_name = 'growthprocess2d'
    sys_call = "./" + exe_name + " " + filename_parameters
    print(sys_call)
    os.system(sys_call)


def homogenize_2d(parameters):
    filename_hom_elas = parameters["name"] + "_tensor.txt"
    print(filename_hom_elas)
    if os.path.exists(filename_hom_elas):
        os.remove(filename_hom_elas)
    octavecall = "octave --no-window-system --no-site-file --no-gui homogenize_2d_image.m"
    syscall = octavecall + " " + parameters["name"] + ".ppm " + filename_hom_elas + " " + str(parameters["young_modulus_solid"]) + " " + str(parameters["poisson_ratio"]) + " " + str(parameters["ignore_void"])
    os.system(syscall)
    assert os.path.exists(filename_hom_elas)
    return numpy.loadtxt(filename_hom_elas)


def compute_porous_material_homogenize_2d_random_seed(num_sample, parameters, random_seed, basename):
    parameters["random_seed"] = str(num_sample + 2 + random_seed)
    parameters["name"] = basename + "_sample" + str(num_sample)
    compute_porous_material(parameters)
    elasticity_tensor = homogenize_2d(parameters)
    os.remove(parameters["name"] + ".ppm")  # save space (remove PPM image)
    return elasticity_tensor


def compute_porous_material_stochastic_homogenize_2d(parameters, parallel=False, max_cores=6):
    elasticity_tensors = []
    num_samples = int(parameters["num_samples_montecarlo"])
    basename = parameters["name"]
    random_seed = int(parameters["random_seed"])
    if not parallel:
        # sequential computing
        for num_sample in range(num_samples):
            elasticity_tensors.append(
                compute_porous_material_homogenize_2d_random_seed(
                    num_sample, parameters, random_seed, basename))
    else:
        # parallel computing
        # recommendation: use export OPENBLAS_NUM_THREADS=1 to disable BLAS threads for improved efficiency
        num_used_cores = min(max(multiprocessing.cpu_count(), 1), max_cores)
        with multiprocessing.Pool(num_used_cores) as pool:
            list_arguments = []
            for num_sample in range(num_samples):
                list_arguments.append((num_sample, parameters, random_seed, basename))
            elasticity_tensors = pool.starmap(compute_porous_material_homogenize_2d_random_seed, list_arguments)
    parameters["name"] = basename
    parameters["random_seed"] = random_seed
    mean_elasticity_tensors = numpy.mean(elasticity_tensors, axis=0)
    print('mean_elasticity_tensors=\n' + str(mean_elasticity_tensors))
    return mean_elasticity_tensors


def parse_compute_porous_material_stochastic_homogenize_2d(parameters):
    # returns the set of elasticity tensors, and results statistics
    elasticity_tensor_set = []
    porosity_set = []
    porosity_regularized_set = []
    num_connected_components_set = []
    num_connected_components_regularized_set = []
    num_sites_set = []
    num_samples = int(parameters["num_samples_montecarlo"])
    basename = parameters["name"]
    random_seed = int(parameters["random_seed"])
    for num_sample in range(num_samples):
        parameters["random_seed"] = str(num_sample + 2 + random_seed)
        parameters["name"] = basename + "_sample" + str(num_sample)
        filename_elasticity_tensor = parameters["name"] + "_tensor.txt"
        filename_results = parameters["name"] + ".res"
        if os.path.isfile(filename_elasticity_tensor) and os.path.isfile(filename_results):
            elasticity_tensor = numpy.loadtxt(parameters["name"] + "_tensor.txt")
            elasticity_tensor_set.append(elasticity_tensor)
            results = read_parameters(filename_results)
            porosity_set.append(float(results["porosity"]))
            porosity_regularized_set.append(float(results["porosity_regularized"]))
            num_connected_components_set.append(float(results["num_connected_components"]))
            num_connected_components_regularized_set.append(float(results["num_connected_components_regularized"]))
            num_sites_set.append(float(results["num_sites"]))
    parameters["name"] = basename
    parameters["random_seed"] = random_seed
    return (elasticity_tensor_set, porosity_set, porosity_regularized_set, num_connected_components_set, num_connected_components_regularized_set, num_sites_set)


def get_closest_isotropic_elasticity_tensor_frobenius(elasticity_tensor):
    c11 = elasticity_tensor[0][0]
    c12 = elasticity_tensor[0][1]
    c22 = elasticity_tensor[1][1]
    c33 = elasticity_tensor[2][2]
    c11_iso = 1.0 / 20.0 * (9.0 * (c11 + c22) + 2.0 * c12 + 4.0 * c33)
    c12_iso = 1.0 / 20.0 * (c11 + c22 + 18.0 * c12 - 4.0 * c33)
    isotropic_elasticity_tensor = numpy.array([
        [c11_iso, c12_iso, 0],
        [c12_iso, c11_iso, 0],
        [0, 0, (c11_iso - c12_iso) / 2.0]])
    # check that the resulting tensor is positive definite
    assert c11_iso > 0
    assert c11_iso - c12_iso > 0
    assert c11_iso * c11_iso - c12_iso * c12_iso > 0
    return isotropic_elasticity_tensor


def get_poissons_ratio(isotropic_elasticity_tensor):
    return isotropic_elasticity_tensor[0][1] / isotropic_elasticity_tensor[0][0]


def get_young_modulus(isotropic_elasticity_tensor):
    return (isotropic_elasticity_tensor[0][0] ** 2 - isotropic_elasticity_tensor[0][1] ** 2) / isotropic_elasticity_tensor[0][0]


def get_deviation_isotropy(elasticity_tensor):
    isotropic_elasticity_tensor = get_closest_isotropic_elasticity_tensor_frobenius(elasticity_tensor)
    return numpy.linalg.norm(elasticity_tensor - isotropic_elasticity_tensor, ord='fro') / numpy.linalg.norm(elasticity_tensor, ord='fro')


def coefficient_variation(x):
    return numpy.std(x) / numpy.mean(x)


def crop_pdf(filename_plot_pdf, filename_plot_cropped_pdf=None):
    if filename_plot_cropped_pdf is None:
        (head, tail) = os.path.split(filename_plot_pdf)
        filename_plot_cropped_pdf = os.path.join(head, "crop_" + tail)
    # crop margins pdf to save space
    os.system("~/.local/bin/pdf-crop-margins -p 0 -o " + filename_plot_cropped_pdf + " " + filename_plot_pdf)  # pip install pdfCropMargins --upgrade --user
    shutil.move(filename_plot_cropped_pdf, filename_plot_pdf)

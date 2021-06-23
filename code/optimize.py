import sys
import numpy
import utils
import matplotlib.pyplot as plt
import skopt  # scikit-optimize
import skopt.plots
from matplotlib import rc
from scipy.optimize import OptimizeResult


# global parameters
parameters = None
lowest_poisson = 1
lowest_poisson_deviation_isotropy = None


def objective_function(input_variables):
    global lowest_poisson
    global lowest_poisson_deviation_isotropy
    print("input variables = " + str(input_variables))
    variables = input_variables
    utils.parse_porous_material_from_variables(parameters, variables)
    elasticity_tensor = utils.compute_porous_material_stochastic_homogenize_2d(parameters)  # stochastic homogenization
    isotropic_elasticity_tensor = utils.get_closest_isotropic_elasticity_tensor_frobenius(elasticity_tensor)
    poisson = utils.get_poissons_ratio(isotropic_elasticity_tensor)
    print ("* Poisson = " + str(poisson))
    deviation_isotropy = utils.get_deviation_isotropy(elasticity_tensor)
    print("* Deviation isotropy = " + str(deviation_isotropy))
    if poisson < lowest_poisson:
        lowest_poisson = poisson
        lowest_poisson_deviation_isotropy = deviation_isotropy
    return poisson


def get_bounds(parameters):
    num_variables = int(parameters["num_radial_spans"])
    min_radial_span = float(parameters["min_radial_span"])
    max_radial_span = float(parameters["max_radial_span"])
    num_variables_distance = num_variables
    if utils.is_lock_max_radial_span(parameters):
        num_variables_distance = num_variables - 1
    lb = min_radial_span * numpy.ones(num_variables_distance)
    ub = max_radial_span * numpy.ones(num_variables_distance)
    if utils.has_parameterized_length_growth(parameters):
        min_length_growth_radial_span = float(parameters["min_length_growth_radial_span"])
        max_length_growth_radial_span = float(parameters["max_length_growth_radial_span"])
        assert max_length_growth_radial_span > min_length_growth_radial_span
        lb = numpy.concatenate((lb, numpy.ones(num_variables) * min_length_growth_radial_span))
        ub = numpy.concatenate((ub, numpy.ones(num_variables) * max_length_growth_radial_span))
    return (lb, ub)


def plot_convergence(results):
    rc('text', usetex=True)
    font = {'family': 'normal', 'size': 32}
    rc('font', **font)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlabel(r"Number of calls", fontsize='large')
    ax.set_ylabel(r"$v$", fontsize='large')
    ax.grid()
    ax.set_ylim(bottom=-1, top=1)
    if isinstance(results, tuple):
        name, results = results
    else:
        name = None
    if isinstance(results, OptimizeResult):
        n_calls = len(results.x_iters)
        mins = [numpy.min(results.func_vals[:i]) for i in range(1, n_calls + 1)]
        ax.plot(range(1, n_calls + 1), mins, c="k", marker="", markersize=2, lw=4, label=name)
        ax.axhline(y=mins[n_calls - 1], linestyle="--", color="r", lw=3)


def optimize_skot_gp_minimize(objective_function, lb, ub, opt_max_eval, opt_n_random_starts):
    # optimization parameters
    random_seed_opt = 42
    dimensions = numpy.stack((lb, ub), axis=-1)
    res = skopt.gp_minimize(
        objective_function,               # the function to minimize
        dimensions,                       # the bounds on each dimension of x
        acq_func="gp_hedge",              # the acquisition function
        n_calls=opt_max_eval,                 # the number of evaluations of f
        n_random_starts=opt_n_random_starts,  # the number of random initialization points
        random_state=random_seed_opt,     # the random seed
        acq_optimizer="lbfgs")
    plot_convergence(res)
    plt.tight_layout()
    filename_plot_pdf = parameters["name"] + "_convergence_plot.pdf"
    filename_plot_cropped_pdf = parameters["name"] + "_convergence_plot_cropped.pdf"
    plt.savefig(filename_plot_pdf)
    utils.crop_pdf(filename_plot_pdf, filename_plot_cropped_pdf)
    return (res.x, res.fun)


if __name__ == "__main__":
    assert len(sys.argv) == 2
    parameters = utils.read_parameters_optimize(sys.argv[1])
    print(parameters)
    (lb, ub) = get_bounds(parameters)
    (optimized_variables, last_optimum_value) = optimize_skot_gp_minimize(
        objective_function, lb, ub, int(parameters["opt_max_eval"]), int(parameters["opt_n_random_starts"]))
    print("optimized_variables = " + str(optimized_variables))
    print("last_optimum_value = " + str(last_optimum_value))
    with open(parameters["name"] + "_last_optimum_value.txt", 'w') as f:
        f.write(str(last_optimum_value))
    print("lowest_poisson = " + str(lowest_poisson))
    if lowest_poisson_deviation_isotropy is not None:
        print("lowest_poisson_deviation_isotropy = " + str(lowest_poisson_deviation_isotropy))
        with open(parameters["name"] + "_last_deviation_isotropy.txt", 'w') as f:
            f.write(str(lowest_poisson_deviation_isotropy))
    print("* Saving last optimum result (image of one random instance)")
    utils.parse_porous_material_from_variables(parameters, optimized_variables, name=parameters["name"] + "_last_optimum")
    utils.compute_porous_material(parameters)

import sys
import numpy
import pandas
import seaborn as sns
from matplotlib.ticker import MaxNLocator

import utils


if __name__ == "__main__":
    # read data
    parameters = utils.read_parameters_single(sys.argv[1])
    (elasticity_tensor_set, porosity_set, porosity_regularized_set, num_connected_components_set, num_connected_components_regularized_set, num_sites_set) = utils.parse_compute_porous_material_stochastic_homogenize_2d(parameters)

    porosity_coefficient_variation = utils.coefficient_variation(porosity_regularized_set)
    print("* Coefficient of variation porosity = " + str(porosity_coefficient_variation))

    with open(parameters["name"] + "_latex.txt", 'w') as f:
        utils.write_latex_command(f, "singleanalysisporositycoeffvar", porosity_coefficient_variation)
        utils.write_latex_command(f, "singleanalysisnumsamplesmc", parameters["num_samples_montecarlo"])

    # plot data
    filename_data_csv = parameters["name"] + "_data.csv"
    file_data_csv = open(filename_data_csv, 'w')
    file_data_csv.write("i,Value,Label\n")
    for num_sample, elasticity_tensor in enumerate(elasticity_tensor_set):
        elasticity_tensor = elasticity_tensor_set[num_sample]
        if num_sample > 0:
            for k in range(num_sample + 1):
                for i in range(3):
                    for j in range(3):
                        if j >= i:
                            file_data_csv.write(str(num_sample) + "," + str(elasticity_tensor_set[k][i][j]) + ",c" + str(i + 1) + str(j + 1) + '\n')
    print("* Deviation isotropy = " + str(utils.get_deviation_isotropy(numpy.mean(elasticity_tensor_set, axis=0))))

    print("* Confidence interval (95%) =")
    C_mean = numpy.mean(elasticity_tensor_set, axis=0)
    C_std = numpy.std(elasticity_tensor_set, axis=0)
    k = float(len(elasticity_tensor_set))
    confidence_interval = 1.96 * numpy.divide(C_std, numpy.sqrt(k))
    print(C_mean)
    print("+-")
    print(confidence_interval)

    file_data_csv.close()
    data_csv = pandas.read_csv(filename_data_csv)

    # plot elasticity tensor data
    sns.set_context("paper")
    sns.set(style="whitegrid", font_scale=2.4, rc={'text.usetex': True})
    sns_plot = sns.relplot(x="i", y="Value", kind='line', hue='Label', aspect=2, legend='brief', data=data_csv)
    sns_plot.fig.autofmt_xdate()
    sns_plot.fig.gca().set_xlabel(r'$i$')
    sns_plot.fig.gca().set_ylabel('')
    sns_plot.fig.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    sns_plot._legend.set_title('')
    # replace labels legend, verify before that they coincide!
    new_labels = ['', r'$c_{11}$', r'$c_{12}$', r'$c_{13}$', r'$c_{22}$', r'$c_{23}$', r'$c_{33}$']
    for t, l in zip(sns_plot._legend.texts, new_labels):
        t.set_text(l)
    filename_plot_pdf = parameters["name"] + "_data.pdf"
    sns_plot.savefig(filename_plot_pdf)
    utils.crop_pdf(filename_plot_pdf)

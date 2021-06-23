import sys
import numpy
import itertools
import pandas
import seaborn as sns

import explore_analysis
import utils


def plot_image_size_field(folder_name, dict_image_size, field_name, field_label, max_field_value, plot_name, boxplot_palette):
    if len(dict_image_size) > 0:
        # plot data
        filename_data_csv = folder_name + "data_" + field_name + ".csv"
        file_data_csv = open(filename_data_csv, 'w')
        file_data_csv.write("image_size," + field_name + "\n")
        for image_size in dict_image_size:
            for field in dict_image_size[image_size]:
                file_data_csv.write(str(image_size) + "," + str(field) + "\n")
        file_data_csv.close()
        data_csv = pandas.read_csv(filename_data_csv)
        sns.set_context("paper")
        sns.set(style="whitegrid", rc={'text.usetex': True})
        ax = sns.boxplot(x="image_size", y=field_name, data=data_csv, whis=numpy.inf, linewidth=0.5, palette=boxplot_palette)
        ax = sns.stripplot(x="image_size", y=field_name, data=data_csv, color=".3", size=1.0)
        ax.figure.autofmt_xdate()
        ax.figure.gca().set_xlabel(r'Size $s$')
        ax.figure.gca().set_ylabel(field_label)
        ax.set_ylim(bottom=0, top=max_field_value)
        ax.set_aspect(1.2 / max_field_value)
        ax.figure.savefig(plot_name)
        ax.clear()
        ax.get_figure().clf()


def dict_values_to_list(d):
    return list(itertools.chain(*list(d.values())))


def plot_dev_iso(folder_name, dev_iso, max_dev_iso):
    filename_plot_pdf = folder_name + "deviation_isotropy.pdf"
    plot_image_size_field(
        folder_name,
        dev_iso,
        field_name="deviation_isotropy",
        field_label=r'$\delta_{iso}$',
        max_field_value=max_dev_iso * max_factor,
        plot_name=filename_plot_pdf,
        boxplot_palette=sns.cubehelix_palette(5, start=2, rot=0.4))
    utils.crop_pdf(filename_plot_pdf)


def plot_por_coeff_var(folder_name, por_coeff_var, max_por_coeff_var):
    filename_plot_pdf = folder_name + "porosity_coefficient_variation.pdf"
    plot_image_size_field(
        folder_name,
        por_coeff_var,
        field_name="porosity_coefficient_variation",
        field_label=r'$c_{v}(p)$',
        max_field_value=max_por_coeff_var * max_factor,
        plot_name=filename_plot_pdf,
        boxplot_palette=sns.cubehelix_palette(5, start=2, rot=0.1))
    utils.crop_pdf(filename_plot_pdf)


if __name__ == "__main__":
    assert(len(sys.argv) == 3)
    (folder_name_1, dev_iso_1, por_coeff_var_1) = explore_analysis.explore_analysis_main(sys.argv[1])
    (folder_name_2, dev_iso_2, por_coeff_var_2) = explore_analysis.explore_analysis_main(sys.argv[2])
    # compute maxium dev_iso and por_coeff_var among 1 and 2 (for the plots)
    max_dev_iso = max(dict_values_to_list(dev_iso_1) + dict_values_to_list(dev_iso_2))
    max_por_coeff_var = max(dict_values_to_list(por_coeff_var_1) + dict_values_to_list(por_coeff_var_2))
    max_factor = 1.2   # relative factor increase for plots
    # plots deviation isotropy
    plot_dev_iso(folder_name_1, dev_iso_1, max_dev_iso)
    plot_dev_iso(folder_name_2, dev_iso_2, max_dev_iso)
    # plots porosity coefficient variation
    plot_por_coeff_var(folder_name_1, por_coeff_var_1, max_por_coeff_var)
    plot_por_coeff_var(folder_name_2, por_coeff_var_2, max_por_coeff_var)

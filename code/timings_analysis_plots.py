import sys
import pandas
import seaborn as sns
import matplotlib.pyplot as plt

import timings_analysis
import utils


def plot_timings_plot(filename_data_csv, plot_name, field_y, field_y_label, y_lim_top=None):
    data_csv = pandas.read_csv(filename_data_csv)
    sns.set(style="whitegrid", rc={'text.usetex': True})
    sns.set_context("paper", font_scale=2.3)
    ax = sns.scatterplot(x="pixels_million", y=field_y, hue='legend_text', data=data_csv, color='gray', edgecolor="none", s=80)
    ax.figure.gca().set_xlabel(r'Number of pixels $n = s \times s$ ($10^6$)')
    ax.figure.gca().set_ylabel(field_y_label)
    ax.figure.gca().set_ylim(bottom=0)
    ax.figure.gca().set_xlim(left=0)
    if y_lim_top is not None:
        ax.figure.gca().set_ylim(top=y_lim_top)
    plt.legend(loc='upper left')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[1:], labels=labels[1:], loc=2)
    ax.figure.savefig(plot_name, bbox_inches='tight')
    ax.clear()
    ax.get_figure().clf()
    utils.crop_pdf(plot_name)


if __name__ == "__main__":
    assert len(sys.argv) > 1
    filename_data_csv = sys.argv[1] + "data_timings.csv"
    file_data_csv = open(filename_data_csv, 'w')
    file_data_csv.write("pixels_million,time_min,legend_text\n")
    with open(sys.argv[1] + "latex.tex", 'w') as file_latex:
        for filename_parameters in sys.argv[2:]:
            print(filename_parameters)
            (folder_name, pixels_million_set, time_growth_process_min_set, range_growth_length) = timings_analysis.timings_analysis_main(filename_parameters, file_latex)
            range_growth_length_legend = None
            if len(range_growth_length) == 1:
                range_growth_length_legend = r'$\psi = ' + str(range_growth_length[0]) + "$"
            elif len(range_growth_length) == 2:
                range_growth_length_legend = r'$\psi \in [' + str(range_growth_length[0]) + "-" + str(range_growth_length[1]) + "]$"
            for i in range(len(pixels_million_set)):
                legend_text = range_growth_length_legend
                file_data_csv.write(str(pixels_million_set[i]) + "," + str(time_growth_process_min_set[i]) + "," + legend_text + "\n")
        file_data_csv.close()

    plot_timings_plot(
        filename_data_csv, sys.argv[1] + "data_timings.pdf",
        field_y="time_min", field_y_label="Time (min)", y_lim_top=12.5)

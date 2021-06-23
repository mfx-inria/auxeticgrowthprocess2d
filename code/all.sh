make clean
make

## figure star-shaped set interpolation
./growthprocess2d parameters/figure_function_interpolation_only_points.txt
./growthprocess2d parameters/figure_function_interpolation.txt
./growthprocess2d parameters/figure_function_interpolation_3fold.txt
./growthprocess2d parameters/figure_function_interpolation_4axis.txt

## figure cell regularization
./growthprocess2d parameters/figure_problem_connectivity.txt

echo "- figures overview growth..." > ../results/log.txt
date >> ../results/log.txt
python3 single.py parameters/figure_overview_euclidean.txt 0
python3 single.py parameters/figure_overview_constant.txt 0
python3 single.py parameters/figure_overview.txt 0

echo "- figure example stochastic homogenization..." >> ../results/log.txt
date >> ../results/log.txt
python3 single.py parameters/figure_example_stochastic_homogenization.txt
python3 single_analysis.py parameters/figure_example_stochastic_homogenization.txt

echo "- exploration three-fold symmetry (RSA)..." >> ../results/log.txt
date >> ../results/log.txt
python3 explore.py parameters/explore_variable_rsa.txt
python3 explore_analysis.py parameters/explore_variable_rsa.txt &> ../results/log_analysis_explore_variable_rsa.txt
echo "- exploration three-fold symmetry (Poisson)..." >> ../results/log.txt
date >> ../results/log.txt
python3 explore.py parameters/explore_variable_random.txt
python3 explore_analysis.py parameters/explore_variable_random.txt &> ../results/log_analysis_explore_variable_random.txt
python3 explore_analysis_plots.py parameters/explore_variable_rsa.txt parameters/explore_variable_random.txt

## minimization Poisson's ratio with different settings
echo "- optimize fixed (Poisson)..." >> ../results/log.txt
date >> ../results/log.txt
python3 optimize.py parameters/optimize_fixed_random.txt
echo "- optimize fixed (RSA)..." >> ../results/log.txt
date >> ../results/log.txt
python3 optimize.py parameters/optimize_fixed_rsa.txt
echo "- optimize variable first (Poisson)..." >> ../results/log.txt
date >> ../results/log.txt
python3 optimize.py parameters/optimize_variable_max_length_growth_first_random.txt
echo "- optimize variable first (RSA)..." >> ../results/log.txt
date >> ../results/log.txt
python3 optimize.py parameters/optimize_variable_max_length_growth_first_rsa.txt
echo "- optimize variable second (Poisson)..." >> ../results/log.txt
date >> ../results/log.txt
python3 optimize.py parameters/optimize_variable_max_length_growth_second_random.txt
echo "- optimize variable second (RSA)..." >> ../results/log.txt
date >> ../results/log.txt
python3 optimize.py parameters/optimize_variable_max_length_growth_second_rsa.txt

# timings Poisson
python3 timings.py ../results/random/optimize_variable_max_length_growth_second/optimize_current_last_optimum.txt
python3 timings.py ../results/random/optimize_variable_max_length_growth_first/optimize_current_last_optimum.txt
python3 timings.py ../results/random/optimize_fixed/optimize_current_last_optimum.txt
python3 timings_analysis_plots.py ../results/timings/poisson_ ../results/random/optimize_fixed/optimize_current_last_optimum.txt ../results/random/optimize_variable_max_length_growth_first/optimize_current_last_optimum.txt ../results/random/optimize_variable_max_length_growth_second/optimize_current_last_optimum.txt

# timings RSA
python3 timings.py ../results/rsa/optimize_variable_max_length_growth_second/optimize_current_last_optimum.txt
python3 timings.py ../results/rsa/optimize_variable_max_length_growth_first/optimize_current_last_optimum.txt
python3 timings.py ../results/rsa/optimize_fixed/optimize_current_last_optimum.txt
python3 timings_analysis_plots.py ../results/timings/rsa_ ../results/rsa/optimize_fixed/optimize_current_last_optimum.txt ../results/rsa/optimize_variable_max_length_growth_first/optimize_current_last_optimum.txt ../results/rsa/optimize_variable_max_length_growth_second/optimize_current_last_optimum.txt

echo "- plot closeup auxetic..." >> ../results/log.txt
date >> ../results/log.txt
python3 plot_closeup_no_video.py parameters/optimize_variable_max_length_growth_second_rsa.txt ../results/auxetic_closeup/auxetic 20000 1000

echo "- All done" >> ../results/log.txt
date >> ../results/log.txt
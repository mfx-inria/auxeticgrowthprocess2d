import sys

import utils

if __name__ == "__main__":
    assert len(sys.argv) >= 2
    parameters = utils.read_parameters_single(sys.argv[1])
    if (len(sys.argv) > 2):
        utils.compute_porous_material(parameters)
    else:
        utils.compute_porous_material_stochastic_homogenize_2d(parameters)

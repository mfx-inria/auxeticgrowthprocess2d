import sys
import PIL
import PIL.ImageFilter
from PIL import Image

import utils


def resample_image(original_image, start, end, max_downsample_image_size, filename):
    newsize = (max_downsample_image_size, max_downsample_image_size)
    radius_gaussian_blur = ((end - start) / (max_downsample_image_size) - 1)
    print(radius_gaussian_blur)
    cropped_image = original_image.filter(PIL.ImageFilter.GaussianBlur(radius=radius_gaussian_blur))
    cropped_image = cropped_image.resize(newsize, resample=PIL.Image.LANCZOS, box=(start, start, end, end))
    cropped_image.save(filename, "png")


if __name__ == "__main__":
    assert len(sys.argv) == 5
    parameters_optimization = utils.read_parameters(sys.argv[1])
    target_name = sys.argv[2]
    target_image_size = int(sys.argv[3])
    max_downsample_image_size = int(sys.argv[4])
    parameters_sample = utils.read_parameters(parameters_optimization["name"] + "_last_optimum.txt")
    parameters_sample["image_size"] = str(target_image_size)
    parameters_sample["name"] = target_name
    parameters_sample["plot_porous_material_sites_png"] = "true"
    parameters_sample["plot_sites_png"] = "false"
    parameters_sample["plot_starshaped_ppm"] = "false"
    parameters_sample["plot_starshaped_pdf"] = "false"
    parameters_sample["save_porous_material_ppm"] = "false"
    parameters_sample["save_results_txt"] = "true"
    parameters_sample["point_process"] = "File"
    parameters_sample["filename_points"] = "data/RSA_s" + str(target_image_size) + "_D" + parameters_sample["rsa_max_dist_rel_pixel"] + ".bin"
    print(parameters_sample["filename_points"])
    utils.compute_porous_material(parameters_sample)
    PIL.Image.MAX_IMAGE_PIXELS = 800000000  # allow to read big images
    original_image = PIL.Image.open(parameters_sample["name"] + "_porous_material.png")

    # paper closeups
    for crop_factor in [1, 2, 5, 10, 15, 20]:
        # centered crop
        crop_size = int(original_image.size[0] / float(crop_factor))
        start = int((original_image.size[0] - crop_size) / 2)
        end = int((original_image.size[0] + crop_size) / 2)
        resample_image(original_image, start, end, max_downsample_image_size, parameters_sample["name"] + "_x" + str(crop_factor) + ".png")

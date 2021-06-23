run("homogenize_2d.m")

function l = lame_first(E, v)
    l = (E * v) / ((1.0 + v) * (1.0 - 2.0 * v));
endfunction

function s = shear_modulus(E, v)
    s = E / (2.0 * (1.0 + v));
endfunction

function ml = modified_lame_first(E, v)
    % Modified Lame's first parameter, to get plane stress properties in 2D.
    lam = lame_first(E, v);
    mu = shear_modulus(E, v);
    ml = (2.0 * mu * lam) / (lam + 2 * mu);
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% usage
% homogenize_2d_image [input_image_filename] [output_tensor_filename] [young_modulus_solid] [poisson_ratio] [ignore_void]
% display(imfinfo(argv(){1}));  % display image information
x = imread(argv(){1});
young_modulus_solid = str2double(argv(){3});
disp("* Young modulus solid = ") , disp(young_modulus_solid);
poissons_ratio = str2double(argv(){4});
disp("* Poisson's ratio = ") , disp(poissons_ratio);
ignore_void = (str2num(argv(){5}) != 0);
if ignore_void
    disp("* Ignoring void material");
    lambda = [0, modified_lame_first(young_modulus_solid, poissons_ratio)];
    mu = [0, shear_modulus(young_modulus_solid, poissons_ratio)];
else
    young_modulus_void = 1e-10;
    lambda = [modified_lame_first(young_modulus_void, poissons_ratio), modified_lame_first(young_modulus_solid, poissons_ratio)];
    mu = [shear_modulus(young_modulus_void, poissons_ratio), shear_modulus(young_modulus_solid, poissons_ratio)];
endif

%% compute homogenized elasticity tensor (using publicly available method of Andreassen and Andreasen 2012)
[nely, nelx] = size(x);
lx = nelx;  % Unit cell length in x-direction.
ly = nely;  % Unit cell length in y-direction.
CH = homogenize_2d_andreassen(lx, ly, lambda, mu, 90, x, ignore_void);
save("-ascii", argv(){2}, "CH");
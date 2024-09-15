# LADAR
SAR LADAR Image Reconstruction

There's a lot of MATLAB files:

generation_smth.m -- generate a data from original image

test/testing_smth.m / RMSD_.m -- code which use the algorithm with some fixed optimal parametres

to_array.m / to_Matrix.m / resize_n.m -- add some help functions

help_omega(_sign).m -- regulariazation matrix to find solution in specific norm

iterative_EM.m -- original approach to reconstruct the image

grad_reg_fun_of_lim_var(_sign).m -- "classical" approach of gradient methods with L1 regularization

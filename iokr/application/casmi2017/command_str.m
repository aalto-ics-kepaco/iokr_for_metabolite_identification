addpath (genpath ('/m/cs/scratch/kepaco/bache1/code/iokr-mp'));
param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), {'debug_param', 'opt_param', 'iokr_param', 'data_param', 'ky_param'});
param.iokr_param.model_representation = 'Chol_decomp_of_C';
param.ky_param.representation  = 'feature';
param.ky_param.type            = 'linear';
param.ky_param.base_kernel     = 'linear';
param.ky_param.param_selection = 'cv';
param.iokr_param.mkl='unimkl';
challenge_param = struct ('ion_mode', 'positive', 'fp_set', 'masked');
input_dir_training='/m/cs/scratch/kepaco/bache1/data/metabolite-identification/CASMI2017/input/';
input_dir_test='/m/cs/scratch/kepaco/bache1/data/metabolite-identification/CASMI2017/input/challenges/';
output_dir='/m/cs/scratch/kepaco/bache1/data/metabolite-identification/CASMI2017/output/iokr/unimkl_linear/';
run_IOKR_CASMI2017 (input_dir_training, input_dir_test, output_dir, param, challenge_param);

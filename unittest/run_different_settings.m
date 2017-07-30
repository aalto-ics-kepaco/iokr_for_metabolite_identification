param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
    {'debug_param', 'opt_param', 'iokr_param', 'data_param', 'ky_param'});
param.ky_param.representation  = 'kernel';
param.ky_param.type            = 'gaussian';
param.ky_param.base_kernel     = 'tanimoto';
param.ky_param.param_selection = 'entropy';
run_IOKR (inputDir, outputDir, cand_transposed__, param)

param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
    {'debug_param', 'opt_param', 'iokr_param', 'data_param', 'ky_param'});
param.ky_param.representation  = 'feature';
param.ky_param.type            = 'linear';
param.ky_param.base_kernel     = 'linear';
param.ky_param.param_selection = 'cv';
run_IOKR (inputDir, outputDir, cand_transposed__, param)

param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
    {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param', 'ky_param'});
param.ky_param.representation  = 'kernel';
param.ky_param.type            = 'gaussian';
param.ky_param.base_kernel     = 'tanimoto';
param.ky_param.param_selection = 'entropy';
run_MP_IOKR (inputDir, outputDir, cand_transposed__, param)

param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
    {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param', 'ky_param'});
param.ky_param.representation  = 'feature';
param.ky_param.type            = 'linear';
param.ky_param.base_kernel     = 'linear';
param.ky_param.param_selection = 'cv';
run_MP_IOKR (inputDir, outputDir, cand_transposed__, param)
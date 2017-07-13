function [approx_error] = aggregate_approximation_error (cv_param, inputKernel, rev_iokr)
    %% Define the input directories
    outputDirMP = ...
        '/m/cs/scratch/kepaco/bache1/data/metabolite-identification/GNPS/results/approx-error/';
    %% Set the parameter
    param = struct ();
    param.mp_iokr_param = struct ('rev_iokr', rev_iokr);
    param.data_param    = struct ('inputKernel', inputKernel, 'cv_param', cv_param); 
    
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, ...
        {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param'});
    
    %% Load the results of MP
    approx_error = zeros (4138, 1);
    settingHash = DataHash (struct (                            ...
        'cv_param',        param.data_param.cv_param,           ...
        'repetition',      param.data_param.repetition,         ...
        'input_kernel',    upper(param.data_param.inputKernel), ...
        'center',          param.mp_iokr_param.center,          ...
        'rev_iokr',        param.mp_iokr_param.rev_iokr));

    resFn = strcat (outputDirMP, '/', settingHash, '.mat'); 
    fprintf ('%s exists %d\n', resFn, exist (resFn, 'file'));

    tmp = load (resFn); 

    approx_error = tmp.result.approx_error;
end % function
function [rank_perc, ranks] = aggregate_results_all_separate (inputKernel, outputDirMP, cv_param, inclExpCand)
    %% Define the input directories
%     outputDirMP = ...
%         '/m/cs/scratch/kepaco/bache1/data/metabolite-identification/GNPS/results/mp-iokr/commit_431c5b22cce549b07374677cbd514ade2c2b63ed/';
    %% Set the parameter
    param = struct ();
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, ...
        {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param'});
    param.mp_iokr_param.rev_iokr = 'separate';

    param.data_param.selection_param = struct ( ...
        'strategy', 'all', 'inclExpCand', inclExpCand);
    
    param.data_param.inputKernel = inputKernel;
    param.data_param.cv_param = cv_param;

    %% Load the results of MP
    rank_perc = zeros (100, 1);
    ranks     = zeros (4138, 1);
    
    settingHash = DataHash (struct (                                ...
        'cv_param',        param.data_param.cv_param,           ...
        'selection_param', param.data_param.selection_param,    ...
        'repetition',      param.data_param.repetition,         ...
        'input_kernel',    upper(param.data_param.inputKernel), ...
        'center',          param.mp_iokr_param.center,          ...
        'rev_iokr',        param.mp_iokr_param.rev_iokr));

    resFn = strcat (outputDirMP, '/', settingHash, '.mat'); 
    fprintf ('%s exists %d\n', resFn, exist (resFn, 'file'));

    tmp = load (resFn); 
    rank_perc(:) = tmp.result.rank_perc(1:100);
    ranks(:) = tmp.result.ranks;
    
    assert (all ( ...
        getRankPerc (tmp.result.ranks, tmp.result.cand_num) == tmp.result.rank_perc));
end % function
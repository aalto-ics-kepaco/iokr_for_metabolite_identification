function [rank_perc, ranks, cand_num] = aggregate_results_maxElement (inputKernel, cv_param, maxNumCand, inclExpCand)
    %% Define the input directories
    outputDirMP = ...
        '/m/cs/scratch/kepaco/bache1/data/metabolite-identification/GNPS/results/mp-iokr/';
    %% Set the parameter
    param = struct ();
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, ...
        {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param'});

    param.data_param.selection_param = struct ( ...
        'strategy', 'maxElement', 'maxNumCand', maxNumCand, 'inclExpCand', inclExpCand);
    
    param.data_param.inputKernel = inputKernel;
    param.data_param.cv_param = cv_param;

    %% Load the results of MP  
    n_rep = 21;
    rank_perc = zeros (100, n_rep);
    ranks     = zeros (4138, n_rep);
    if (nargout > 2)
        cand_num = zeros (4138, n_rep);
    end % if
    
    for rep = 1:n_rep      
        param_rep = param;
        param_rep.data_param.repetition = rep;

        settingHash = DataHash (struct (                                ...
            'cv_param',        param_rep.data_param.cv_param,           ...
            'selection_param', param_rep.data_param.selection_param,    ...
            'repetition',      param_rep.data_param.repetition,         ...
            'input_kernel',    upper(param_rep.data_param.inputKernel), ...
            'center',          param_rep.mp_iokr_param.center,          ...
            'rev_iokr',        param_rep.mp_iokr_param.rev_iokr));

        resFn = strcat (outputDirMP, '/', settingHash, '.mat'); 
%         fprintf ('%s exists %d\n', resFn, exist (resFn, 'file'));
        
        tmp = load (resFn); 
        rank_perc(:, rep) = tmp.result.rank_perc(1:100);
        ranks(:, rep) = tmp.result.ranks;
        
        if (nargout > 2)
            cand_num(:, rep) = tmp.result.cand_num_sel;
        end % if    
        
        assert (all ( ...
            getRankPerc (tmp.result.ranks, tmp.result.cand_num) == tmp.result.rank_perc));
    end % for
end % function
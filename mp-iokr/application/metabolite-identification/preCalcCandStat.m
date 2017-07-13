function preCalcCandStat (param, useAllExamplesForTraining)
    %% Set parameter
    if (nargin < 1)
        param = struct ();
    end % if
    if (nargin < 2)
        useAllExamplesForTraining = false;
    end % if
        
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, ...
        {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param'});
    %% Load data    
    inputDir = '/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/';
    inOutDir = strcat (inputDir, '/pre_calculated_stats/');
    if (param.debug_param.isDebugMode)
        inOutDir = strcat (inOutDir, '/debug/');
    end % if
    
    % Fingerprints
    Y = load (strcat (inputDir, '/fingerprints/fp.mat'));
    Y = full (Y.Y);
    
    rng (param.debug_param.randomSeed);
    
    if (param.debug_param.isDebugMode)
        fprintf ('Show me a random-number: %f\n', rand(1));
        
        % The debug-mode is mainly used to have some dry run and check all
        % the functionality. Nothing should explode during that run, which
        % gives evidence that in the use-case we also do not encounter
        % problems. Therefore lets make the following simplifications:
        param.opt_param = struct (   ...
            'val_gamma',   [0.5, 1], ...
            'val_lambda',  [0.5, 1], ...
            'nOuterFolds', 10,       ...
            'nInnerFolds', 2);        
        
        n = size (Y, 2);
        param.debug_param.debug_set = false (n, 1);        
        param.debug_param.debug_set(randsample (n, param.debug_param.n_debug_set)) = true;
        
        Y = Y(:, param.debug_param.debug_set);
    end % if

    % Identifier (inchis) for the training examples
    inchis = readtext (strcat (inputDir, '/inchi.txt'));
    if (param.debug_param.isDebugMode)
        inchis = inchis(param.debug_param.debug_set);
    end % if
    
    % Pre-defined cross-validation folds to prevent to train using
    % molecular structures which are contained in the test
    ind_fold = load (strcat (inputDir, '/cv_ind.txt'));
    if (param.debug_param.isDebugMode)
        ind_fold = ind_fold(param.debug_param.debug_set);
    end % if

    % Candidate sets 
    mf_corres = load (strcat (inputDir, '/candidates/matching_mf_train.txt'));
    if (param.debug_param.isDebugMode)
        mf_corres = mf_corres(param.debug_param.debug_set);
    end % if
    if (~ isempty (param.debug_param.cand))
        warning ('Candidate taken from PARAM.DEBUG_PARAM.CAND');
        
        Y_C = CandidateSets (DataHandle (param.debug_param.cand), mf_corres);      
    else
        cand = load (strcat (inputDir, '/candidates/GNPS_cand_as_struct_transp.mat'));
        Y_C = CandidateSets (DataHandle (cand.cand_struct), mf_corres);
        
        clear cand;
    end 
    
    %% Pre-calculate statistics
    % Cross-validation settings
    if (useAllExamplesForTraining)
        cv_param = struct ('nObservations', numel (ind_fold),               ...
                           'outer', struct ('type', 'random', 'nFolds', 1), ...
                           'inner', struct ('nFolds', param.opt_param.nInnerFolds));   
    else
        cv_param = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold), ...
                           'inner', struct ('nFolds', param.opt_param.nInnerFolds));                       
    end % if
    param.data_param.cv_param = cv_param;
    
    tic;
    getPreCalcCandStat_feat (Y, Y_C, inchis, param, inOutDir);
    fprintf ('Loading / pre-calculating of the candidate statistics took %.3fs\n', toc);
end % if
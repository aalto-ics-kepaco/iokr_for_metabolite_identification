function IOKR_MP_reverse_evaluation (inputDir, outputDir, param)
%% IOKR_MP_EVALUATION evaluation of the metabolite-identification using IOKR-MP
%    ALGORITHM:
%       = Load data associated with the metabolites = 
%           (1) load input-kernels (--> list of kernels for MKL)
%           (2) load fingerprints
%           (3) load candidate-sets
%
%       = Precalculate / load scenarios and statistics =
%           SEE ALGO

% create a candidate selection
% create cross-validation indices (outer and inner)
% pre-calculate the covariance matricies (using cv indices)

% LOOP outer cross-validation
%   

% ALGO
%   - provide selection paratemer
%   - provide cv-partition for the outer fold
%
%   - IF    pre-calculated data should be used --> load data from drive
%   - ELSE  calcuate: cv-partition, selection and statistic

    %% Check the input arguments and set defaults
    if (nargin < 2)
        error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', ...
            'Not enough input arguments.');
    end % if
    if (nargin < 3)
        param = struct ();
    end % if   
    if (~ exist (inputDir, 'dir'))
        error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', '%s: No such directory.', ...
            inputDir);
    end % if
    if (~ exist (outputDir, 'dir'))
        error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', '%s: No such directory.', ...
            outputDir);
    end % if
    
    % Set the defaults values for the parameter in PARAM if needed.
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, ...
        {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param', 'ky_param'});    
    
    assert (~ param.data_param.usePreCalcStat, ...
        'Not yet implemented for the new train and test methods.');
    
    %% Load data 
    % ... input-kernels for the training examples
    if (param.debug_param.isDebugMode)
        % Lets reduce the size of the UNIMKL a bit. This is useful if we
        % want to debug the 'separate' kernel combination.
        param.data_param.availInputKernels = {'PPKR', 'NSF', 'CEC', 'CPJ'};
    end % if
    [KX_list, param] = loadInputKernelsIntoList (strcat (inputDir, '/input_kernels/'), param);
    if (isempty (KX_list))
        error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', ...
            'No kernel loaded.');
    end % if
    
    % Set seed in order to be able to reproduce the settings.
    % By default this is set to 'shuffle'. This means, that the current
    % time is used as seed. However, if we run the script on triton the
    % jobs might start at the same time. We therefore need to be able to
    % set the seed manualy. For example by using the slurm-job-id. This can
    % be done by the calling sbatch-file.
    rng (param.debug_param.randomSeed);
    
    if (param.debug_param.isDebugMode)
        fprintf ('Show me a random-number: %f\n', rand(1));
        
        % Modify the output-dir & creat it if needed
        outputDir = strcat (outputDir, '/debug/');
        if (~ exist (outputDir, 'dir'))
            if (~ mkdir (outputDir))
                error ('IOKR_MP_reverse_feat_evaluation:RuntimeError', ...
                    'Could not create debug-directory: %s', outputDir);
            end % if
        end % if
        
        % The debug-mode is mainly used to have some dry run and check all
        % the functionality. Nothing should explode during that run, which
        % gives evidence that in the use-case we also do not encounter
        % problems. Therefore lets make the following simplifications:
        param.opt_param = struct (   ...
            'val_gamma',   [0.5, 1], ...
            'val_lambda',  [0.5, 1], ...
            'nOuterFolds', 10,       ...
            'nInnerFolds', 2);        
        
        n = size (KX_list{1}, 1);
        param.debug_param.debug_set = false (n, 1);        
        param.debug_param.debug_set(randsample (n, param.debug_param.n_debug_set)) = true;
        KX_list = cellfun(@(x) x(param.debug_param.debug_set, param.debug_param.debug_set), KX_list, 'UniformOutput', false);
    end % if
    
    % ... fingerprints for the training examples
    Y = load (strcat (inputDir, '/fingerprints/fp.mat'));
    Y = full (Y.Y);
    if (param.debug_param.isDebugMode)
        Y = Y(:, param.debug_param.debug_set);
    end % if
    
    % ... identifier (inchis) for the training examples
    inchis = readtext (strcat (inputDir, '/inchi.txt'));
    if (param.debug_param.isDebugMode)
        inchis = inchis(param.debug_param.debug_set);
    end % if
    
    % ... pre-defined cross-validation folds to prevent to train using
    % molecular structures which are contained in the test
    ind_fold = load (strcat (inputDir, '/cv_ind.txt'));
    if (param.debug_param.isDebugMode)
        ind_fold = ind_fold(param.debug_param.debug_set);
    end % if
    
    % ... candidate sets 
    mf_corres = load (strcat (inputDir, '/candidates/matching_mf_train.txt'));
    if (param.debug_param.isDebugMode)
        mf_corres = mf_corres(param.debug_param.debug_set);
    end % if
    tic;
    if (~ isempty (param.debug_param.cand))
        warning ('Candidate taken from PARAM.DEBUG_PARAM.CAND');
        
        Y_C = CandidateSets (DataHandle (param.debug_param.cand), mf_corres);
        
        param.debug_param.cand = []; 
    else
        cand = load (strcat (inputDir, '/candidates/GNPS_cand_as_struct_transp.mat'));
        Y_C = CandidateSets (DataHandle (cand.cand_struct), mf_corres);
        
        clear cand;
    end 
    toc;
    
    ind_eval = load (strcat (inputDir, '/candidates/ind_eval.txt'));
    eval_set = false (numel (mf_corres), 1);
    eval_set(ind_eval) = true;
    if (param.debug_param.isDebugMode)
        eval_set = eval_set(param.debug_param.debug_set);
    end % if
    
    assert (Y_C.getNumberOfExamples() == size (Y, 2), ...
        'The number of examples in the associated with the candidate sets is different from the number of example fingerprint vectors.');    
    
    %% Load / Store pre-calculated statistics
    if (param.data_param.usePreCalcStat)
        cv_param = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold), ...
                           'inner', struct ('nFolds', param.opt_param.nInnerFolds));                       
        param.data_param.cv_param = cv_param;
        
        % NOTE: The selec_ property of Y_C will be modified according to
        %       the selection defined by PARAM.DATA_PARAM.SELECTION_PARAM.
        tic;
        matObj = getPreCalcCandStat_feat (Y, Y_C, inchis, param, strcat (inputDir, '/pre_calculated_stats/'));
        fprintf ('Loading / pre-calculating of the candidate statistics took %.3fs\n', toc);
        
        param.data_param.cv         = matObj.cv;
        param.mp_iokr_param.center  = matObj.center;
        param.data_param.repetition = matObj.repetition;
        % matObj also contains the statistics
        param.data_param.matObj     = matObj;
    else       
        % Select a subset of candidates. By using the default values ALL the 
        % candidates are selected.
        selec = getCandidateSelection (Y_C, inchis, param.data_param.selection_param);      
        Y_C.setSelectionsOfCandidateSets (selec);
        
        cv_param = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold));
        param.data_param.cv_param = cv_param;
        param.data_param.cv       = getCVIndices (param.data_param.cv_param);
    end % if
    
    %% Evaluate the performance using 10-fold cv
    ranks = NaN (Y_C.getNumberOfExamples(), 1);

    for foldIdx = 1:param.data_param.cv.outer.NumTestSets
        fprintf ('Outer fold: %d/%d\n', foldIdx, param.data_param.cv.outer.NumTestSets);
              
        data_param_fold = struct();
        data_param_fold.train_set      = training_my (param.data_param.cv.outer, foldIdx);
        data_param_fold.test_set       = test_my (param.data_param.cv.outer, foldIdx);
        data_param_fold.usePreCalcStat = param.data_param.usePreCalcStat;
        
        if (data_param_fold.usePreCalcStat)
            % By loading the covariance and mean vectors at this place, we
            % save some memory, as only the data is loaded which is needed.
            data_param_fold.cv       = param.data_param.cv.inner{foldIdx};
            data_param_fold.stats    = param.data_param.matObj.stats(1, foldIdx);
            data_param_fold.stats_cv = param.data_param.matObj.stats_cv(:, foldIdx);
        end % if
        
        % Calculate the scores for each candidate corresponding the current
        % test-examples
        % Training
        KX_list_train = cellfun(@(x) x(data_param_fold.train_set, data_param_fold.train_set), ...
            KX_list, 'UniformOutput', false);
        Y_train       = Y(:, data_param_fold.train_set);
        Y_C_train     = Y_C.getSubset (data_param_fold.train_set);
        trained_model = Train_MPIOKR (KX_list_train, Y_train, Y_C_train, ...
            param.ky_param, param.mp_iokr_param, param.opt_param, param.debug_param);
        
        % Scoring
        KX_list_train_test = cellfun(@(x) x(data_param_fold.train_set, data_param_fold.test_set), ...
            KX_list, 'UniformOutput', false);
        KX_list_test       = cellfun(@(x) x(data_param_fold.test_set, data_param_fold.test_set),  ...
            KX_list, 'UniformOutput', false);
        Y_C_test           = Y_C.getSubset (data_param_fold.test_set);
        scores_test = Test_MPIOKR (KX_list_train_test, KX_list_test, trained_model, ...
            Y_train, Y_C_train, Y_C_test, param.mp_iokr_param, param.mp_iokr_param.center, param.debug_param);
        assert (numel (scores_test) == numel (find (data_param_fold.test_set)), ...
            'There must be a score for each test-examples.');
        
        % Computation of the ranks of the test examples
        inchis_test   = inchis(data_param_fold.test_set);
        eval_set_test = eval_set(data_param_fold.test_set);
        
        ranks_test = getRanksBasedOnScores (Y_C_test, inchis_test, scores_test, eval_set_test);
        ranks(data_param_fold.test_set) = ranks_test;
        
        clear data_param_fold;
    end % for

    assert (all (isnan (ranks) == (~ eval_set)), ...
        'The examples without propper candidate set and examples without rank (due to the absence if a candidate set) must be equal.');
    
    %% Calculate rank percentages
    candNum = arrayfun (@(idx) Y_C.getCandidateSet (idx, 0, 'num'), 1:Y_C.getNumberOfExamples());
    rankPerc = getRankPerc (ranks, candNum);
    
    %% Store results
    % FIXME: This result should be stored using the hash value for the
    % corresponding setting. 
    candNumSel = arrayfun (@(idx) Y_C.getCandidateSet (idx, 1, 'num'), 1:Y_C.getNumberOfExamples());
    
    result = struct ('ranks', ranks, 'rank_perc', rankPerc, 'cand_num', candNum, ...
        'cand_num_sel', candNumSel, 'param', param); 
    
    settingHash = DataHash (struct (                            ...
        'cv_param',        param.data_param.cv_param,           ...
        'selection_param', param.data_param.selection_param,    ...
        'repetition',      param.data_param.repetition,         ...
        'input_kernel',    upper(param.data_param.inputKernel), ...
        'center',          param.mp_iokr_param.center,          ...
        'rev_iokr',        param.mp_iokr_param.rev_iokr,        ...
        'ky_param',        param.ky_param));
    save (strcat (outputDir, '/', settingHash, '.mat'), 'result', '-v7.3');
    
    disp (rankPerc([1, 5, 10, 20]));
    disp ('! Ready !');
end % function 
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

function IOKR_MP_reverse_feat_evaluation (inputDir, outputDir, param)
    %% Check the input arguments and set defaults
    if (nargin < 2)
        error ('MP_IOKR_Defaults:setDefaultsIfNeeded:InvalidInput', ...
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
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param'});
    
    if (param.debug_param.isDebugMode)
        rng (param.debug_param.randomSeed);
    end % if
    %% Load data 
    % ... input-kernels for the training examples
    [KX_list, param] = loadInputKernelsIntoList (inputDir, param);
    
    % ... fingerprints for the training examples
    Y = load (strcat (inputDir, '/fingerprints/fp.mat'));
    Y = Y.Y;
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
    if (param.debug_param.isDebugMode)
        cand = param.debug_param.cand_struct; 
    else
        cand = load (strcat (inputDir, '/candidates/GNPS_cand_as_struct_transp.mat'));
        cand = cand.cand_struct;
    end 
    toc;
    Y_C = CandidateSets (DataHandle (cand), mf_corres);
    
    if (param.debug_param.isDebugMode)
        ind_eval = load (strcat (inputDir, '/candidates/ind_eval.txt'));
        eval_set = false (4138, 1);
        eval_set(ind_eval) = true;
        eval_set = eval_set(param.debug_param.debug_set);
    else
        ind_eval = load (strcat (inputDir, '/candidates/ind_eval.txt'));
        eval_set = false (numel (mf_corres), 1);
        eval_set(ind_eval) = true;
    end % if
    
    clear cand;
      
    assert (Y_C.getNumberOfExamples() == size (Y, 2), ...
        'The number of examples in the associated with the candidate sets is different from the number of example fingerprint vectors.');
    
    %% Load / Store pre-calculated statistics
    if (param.data_param.usePreCalcStat)
        cv_param = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold), ...
                           'inner', struct ('nFolds', param.opt_param.nInnerFolds));                       
        param.opt_param.cv_param = cv_param;
        
        % NOTE: The selec_ property of Y_C will be modified according to
        %       the selection defined by PARAM.DATA_PARAM.SELECTION_PARAM.
        tic;
        matObj = getPreCalcCandStat_feat (Y, Y_C, inchis, param, outputDir);
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
        param.data_param.cv = getCVIndices (cv_param);
    end % if
   
    
    %% Evaluate the performance using 10-fold cv
    ranks = NaN (Y_C.getNumberOfExamples(), 1);
    
    for i = 1:param.data_param.cv.outer.NumTestSets
        fprintf ('Outer fold: %d/%d\n', i, param.data_param.cv.outer.NumTestSets);
              
        data_param_fold = struct();
        data_param_fold.train_set = training_my (param.data_param.cv.outer, i);
        data_param_fold.test_set  = test_my (param.data_param.cv.outer, i);
        data_param_fold.usePreCalcStat = param.data_param.usePreCalcStat;
        
        if (data_param_fold.usePreCalcStat)
            % By loading the covariance and mean vectors at this place, we
            % save some memory, as only the data is loaded which is needed.
            data_param_fold.cv       = param.data_param.cv.inner{i};
            data_param_fold.stats    = param.data_param.matObj.stats(i, 1);
            data_param_fold.stats_cv = param.data_param.matObj.stats_cv(i, :);
        end % if
        
        % Calculate the scores for each candidate corresponding the current
        % test-examples
        tic; 
        
        scores = MP_IOKR_reverse_feat (KX_list, Y(:, data_param_fold.train_set), Y_C, ...
            param.opt_param, param.mp_iokr_param, data_param_fold);
        
        fprintf ('Calculating the scores took %.3fs\n', toc); 
        
        % Computation of the ranks of the test examples
        Y_C_test      = Y_C.getSubset (data_param_fold.test_set);
        inchis_test   = inchis(data_param_fold.test_set);
        ranks_test    = NaN (Y_C_test.getNumberOfExamples(), 1);
        eval_set_test = eval_set (data_param_fold.test_set);
        
        for j = 1:Y_C_test.getNumberOfExamples()
            % In the case of the metabolite identification we should only
            % evaluate the examples from the eval-set. This ensures
            % comparable results with the recent publication of Celine
            % [Brouard2016].
            if (~ eval_set_test(i))
                continue;
            end % if
            
            % Get the inchis of all the candidate in the set for test
            % example j
            inchis_c = Y_C_test.getCandidateSet (j, 0, 'id');
            
            % FIXME: This is just for testing.
            mf_corres_test = mf_corres(data_param_fold.test_set);
            num = Y_C_test.getCandidateSet (j, 0, 'num');
            assert (isnan (num) == isnan (mf_corres_test(j)));
            
            % There might be scores which are NaN, as not each test example
            % has a candidate set. 
            if (isnan (scores{j})) ; continue ;  end % if
            
            [~ , IX] = sort (scores{j}, 'descend');
            
            ind = find (strcmp (inchis_c(IX), inchis_test(j)));
            if (isempty (ind))
                warning ('The inchi of the example itself is not in the candidate set.');               
                continue; 
            end % if
            
            ranks_test(j) = find (strcmp (inchis_c(IX), inchis_test{j}));
        end % for
        
        ranks(data_param_fold.test_set) = ranks_test;
        
        clear data_param_fold;
    end % for
%     
    assert (all (isnan (ranks) == (~ eval_set)), ...
        'The examples without propper candidate set and examples without rank (due to the absence if a candidate set) must be equal.');
    
    candNum = arrayfun (@(idx) Y_C.getCandidateSet (idx, 0, 'num'), 1:Y_C.getNumberOfExamples());
    candNum(isnan (ranks)) = NaN;
    
    rankPerc = getRankPerc (ranks, max (candNum));
    
    % FIXME: This result should be stored using the hash value for the
    % corresponding setting. 
    result = struct ('ranks', ranks, 'rank_perc', rankPerc, 'cand_num', candNum);
    save (strcat (outputDir, '/first_results.mat'), 'result', '-v7.3');
    
    disp (ranks);
    disp (rankPerc);
end % function 

function [KX_list, param] = loadInputKernelsIntoList (inputDir, param)
%% LOADINPUTKERNELLIST
    switch (upper (param.data_param.inputKernel))
        case upper (param.data_param.availInputKernels)
            disp (['Evaluation using a single kernel: ', param.data_param.inputKernel]);

            KX_list = { loadKernel(strcat (inputDir, '/input_kernels/', upper (param.data_param.inputKernel), '.mat')) };

            % If a single kernel is used, than we force the "mkl" option to
            % be 'unimkl'. The weight for the kernel will be 1.
            param.mp_iokr_param.mkl = 'unimkl';
        case upper ({'unimkl', 'alignf'})
            disp (['Evaluation using multiple kernel learning: ', param.data_param.inputKernel]);

            KX_list = cellfun (@(kernelName) loadKernel (strcat (inputDir, '/input_kernels/', kernelName, '.mat')), ...
                upper (param.data_param.availInputKernels), 'UniformOutput', false);

            param.mp_iokr_param.mkl = lower (param.data_param.inputKernel);
        otherwise
            error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', ...
                '%s: Not a valid input kernel. See MP_IOKR_DEFAULTS for a list of available kernels.', ...
                param.data_param.inputKernel);
    end % switch
    
    if (param.debug_param.isDebugMode)
        n = size (KX_list{1}, 1);
        param.debug_param.debug_set = false (n, 1);        
        param.debug_param.debug_set(randsample (n, param.debug_param.n_debug_set)) = true;
        KX_list = cellfun(@(x) x(param.debug_param.debug_set, param.debug_param.debug_set), KX_list, 'UniformOutput', false);
    end % if
end % function

function rankPerc = getRankPerc (ranks, maxCandNum)
%% GETRANKPERC Calculate percentage for the different rankings.
%    INPUT:
%       ranks         ranking for each molecule
%       maxCandNum    maximum number of candidates for a molecule
%
%    OUTPUT:
%       rankPerc     percantages for rank 1 - 100
    nValidRanks = sum (~ isnan (ranks));

    nel = hist (ranks, 1:max(maxCandNum));
    rankPerc = cumsum(nel)';
    rankPerc = rankPerc / nValidRanks * 100;
%     rankPerc = rankPerc(1:100); 
end % function

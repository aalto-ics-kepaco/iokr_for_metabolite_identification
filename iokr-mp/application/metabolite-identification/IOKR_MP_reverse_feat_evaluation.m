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
    cand = load (strcat (inputDir, '/candidates/GNPS_cand_as_struct_transp.mat'));
    toc;
    cand = cand.cand_struct;
    Y_C = CandidateSets (DataHandle (cand), mf_corres);
    
    clear cand mf_corres;
      
    assert (Y_C.getNumberOfExamples() == size (Y, 2), ...
        'The number of examples in the associated with the candidate sets is different from the number of example fingerprint vectors.');
    
    %% Load / Store pre-calculated statistics
    if (param.data_param.usePreCalcStat)
        % If precaclulated statistics are used, the following things have
        % to be loaded:
        %   1.1) cross-validation data splitting
        %   1.2) candidate selection 
        %   (depending on 1.1 & 1.2) --> covariance matrices & mean vectors
        % Just to be sure: All these three things will be in ONE file. 
        
        % Calculate a hash value using the cross-validation and candidate
        % selection parameter. 
        cv_param = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold), ...
                           'inner', struct ('nFolds', param.opt_param.nInnerFolds));
        
        % Use hash to find file matching the parameters
        statHash = DataHash (struct ('cv_param',        cvParam,                          ...
                                     'selection_param', param.data_param.selection_param, ...
                                     'repetition',      param.data_param.repetition,      ...
                                     'center',          param.mp_iokr_param.center));
        statFn = strcat (inputDir, '/pre_calculated_stats/', statHash, '.mat');
        
        if (exist (statFn, 'file'))
            disp ('Load pre-caculated statistics from file: ', statFn);
            
            param.data_param.statsMatObj = matfile (statFn);
            
            % Use pre-calculated statistics to:
            % ... fix the outer and inner folds 
            param.data_param.cv = param.data_param.statsMatObj.cv;
            assert (isequal (param.data_param.cv, getCVIndices (cv_param)), ... 
                'The expected and the loaded CV_PARAM do not match.');
            % ... select candidates 
            selec = param.data_param.statsMatObj.selec;
            % ... choose kernel centralization (true/false)
            if (param.mp_iokr_param.center ~= param.data_param.statsMatObj.center)
                warning ('Kernel centering is always according to the pre-calculated statistics: center = %d.', ...
                    param.data_param.statsMatObj.center);
            end % if
            param.mp_iokr_param.center = param.data_param.statsMatObj.center;
            % ... provide the covariance matrices and mean vectors 
            % Will be loaded later ...
        else
            warning ('No pre-calculated statistics available. %s: No such file.', ...
                statFn);
            disp ('Pre-calculate statistics');
            
            % Select a subset of candidates.
            selec = getCandidateSelection (Y_C, inchis, param.data_param.selection_param); % Can carry randomness.
            param.data_param.cv = getCVIndices (cv_param);                                 % Carries randomness.
            
            try 
                % Create a new statistic file
                matObj = matfile (statFn, 'Writable', true);
                
                matObj.cv         = param.data_param.cv;
                matObj.selec      = selec;
                matObj.repetition = param.data_param.repetition;
                matObj.ker_center = param.mp_iokr_param.center;
                matObj.stats      = struct (1, cv.outer.NumTestSets);
                matObj.stats_cv   = struct (cv.inner.NumTestSets, cv.outer.NumTestSets);
                
                for i = 1:cv.outer.NumTestSets
                    train_set = training_my (param.data_param.cv.outer, i);
                    Y_train = Y(:, train_set);
                    mean_Y_train = mean (Y_train, 2);
                    Y_C_train = Y_C.getSubset (train_set);
                    
                    [Mean_Psi_C_train, Cov_Psi_C_train] = Compute_cov_mean_feat ( ...
                        Y_C_train, mean_Y_train, param.mp_iokr_param.center);     
                    
                    matObj.stats(1, i) = { struct('Mean_Psi_C_train', Mean_Psi_C_train, ...
                                                  'Cov_Psi_C_train',  Cov_Psi_C_train) }; 
                                              
                    clear Mean_Psi_C_train Cov_Psi_C_train;
                    
                    for j = 1:cv.inner.NumTestSets
                        train_set_cv = training_my (param.data_param.cv.inner{i}, j);
                        Y_train_cv = Y_train(:, train_set_cv);
                        mean_Y_train_cv = mean (Y_train_cv, 2);
                        Y_C_train_cv = Y_C_train.getSubset (train_set_cv);
                        
                        [Mean_Psi_C_train_cv, Cov_Psi_C_train_cv] = Compute_cov_mean_feat ( ...
                            Y_C_train_cv, mean_Y_train_cv, param.mp_iokr_param.center);     
                        
                        matObj.stats_cv(j, i) = { struct('Mean_Psi_C_train_cv', Mean_Psi_C_train_cv, ...
                                                         'Cov_Psi_C_train_cv',  Cov_Psi_C_train_cv) }; 
                                                     
                        clear Mean_Psi_C_train_cv Cov_Psi_C_train_cv train_set_cv Y_train_cv mean_Y_train_cv Y_C_train_cv;
                    end % for
                end % for
                
                clear train_set Y_train mean_Y_train Y_C_train; 
            catch me
                disp (me.message);
                
                % If anything goes wrong the statistic file should be
                % removed
                delete (statFn);
            end % try           
        end % if    
    else       
        % Select a subset of candidates. By using the default values ALL the 
        % candidates are selected.
        selec = getCandidateSelection (Y_C, inchis, param.data_param.selection_param);      
        
        cv_param = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold));
        param.data_param.cv = getCVIndices (cv_param);
    end % if
    
    Y_C.setSelectionsOfCandidateSets (selec);
    
    %% Evaluate the performance using 10-fold cv
    ranks = NaN (Y_C.getNumberOfExamples(), 1);
    
    nOuterFolds = param.data_param.cv.outer.NumTestSets;
    for i = 1:nOuterFolds
        data_param_fold = struct();
        data_param_fold.train_set = training_my (param.data_param.cv.outer, i);
        data_param_fold.test_set  = test_my (param.data_param.cv.outer, i);
        data_param_fold.usePreCalcStat = param.data_param.usePreCalcStat;
        
        if (data_param_fold.usePreCalcStat)
            % By loading the covariance and mean vectors at this place, we
            % save some memory, as only the data is loaded which is needed.
            data_param_fold.cv       = param.data_param.cv.inner;
            stats                    = param.data_param.statsMatObj.stats(i, 1);
            data_param_fold.stats    = stats{1};
            stats_cv                 = param.data_param.statsMatObj.stats_cv(i, :);
            data_param_fold.stats_cv = stats_cv{1};
        end % if
          
        scores = MP_IOKR_reverse_feat (KX_list, Y(:, data_param_fold.train_set), Y_C, ...
            param.opt_param, param.mp_iokr_param, data_param_fold);
        
        % Computation of the ranks of the test examples
        Y_C_test    = Y_C.getSubset (data_param_fold.test_set);
        inchis_test = inchis(data_param_fold.test_set);
        ranks_test  = NaN (Y_C_test.getNumberOfExamples(), 1);
        
        for j = 1:Y_C_test.getNumberOfExamples()
            % Get the inchis of all the candidate in the set for test
            % example j
            inchis_c = Y_C_test.getCandidateSet (j, 0, 'id');
            if (~ isnan (inchis_c)) ; continue ; end % if
            
            [~ , IX] = sort (scores{j}, 'descend');
            
            ranks_test(j) = find (strcmp (inchis_c(IX), inchis_test{j}));
        end % for
        
        ranks(data_param_fold.test_set) = ranks_test;
        
        clear data_param_fold;
    end % for
    
    assert (all (isnan (ranks) == isnan (mf_corres)), ...
        'The examples without candidate set and examples without rank (due to the absence if a candidate set) must be equal.');
    
    maxCandNum = arrayfun (@(idx) Y_C.getCandidateSet (idx, 0, 'num'), 1:Y_C.getNumberOfExamples());
    
    rankPerc = getRankPerc (ranks, maxCandNum);
    
    disp (rank);
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
    rankPerc = rankPerc(1:100); 
end % function

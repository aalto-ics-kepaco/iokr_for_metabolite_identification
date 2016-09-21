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
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, {'opt_param', 'mp_iokr_param', 'data_param'});
    
    %% Load data 
    % ... input-kernels for the training examples
    [KX_list, param] = loadInputKernelList (inputDir, param);
    
    % ... fingerprints for the training examples
    Y = load (strcat (inputDir, '/fingerprints/fp.mat'));
    Y = Y.Y;
    
    % ... identifier (inchis) for the training examples
    inchis = readtext (strcat (inputDir, '/inchi.txt'));
    
    % ... candidate sets 
    mf_corres = load (strcat (inputDir, '/candidates/matching_mf_train.txt'));
    cand = load (strcat (inputDir, '/candidates/GNPS_cand.mat'));
    cand = cand.cand;
    Y_C = CandidateSets (DataHandle (cand), mf_corres);
    
    assert (Y_C.getNumberOfExamples() == size (Y, 2), ...
        'The number of examples in the associated with the candidate sets is different from the number of example fingerprint vectors.');
    
    clear cand mf_corres;
    
    %% Evaluate the performance using 10-fold cv
    ranks = NaN (Y_C.getNumberOfExamples(), 1);
    
    if (param.data_param.usePreCalcStat)
       
    else
        ind_fold = load (strcat (inputDir, '/cv_ind.txt'));
        cvParam = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold));
        cv = getCVIndices (cvParam);
    end % if
    
    for i = 1:cv.outer.NumTestSets
        param_cv = param; 
        param_cv.train_set = training_my (cv, i);
        param_cv.test_set  = test_my (cv, i);
        
        scores = MP_IOKR_reverse_feat (KX_list, Y(:, param_cv.train_set), Y_C, ...
            param.opt_param, param.mp_iokr_param, param.data_param);
        
        % Computation of the ranks of the test examples
        Y_C_test = Y_C.getSubset (param_cv.test_set);
        inchis_test = inchis(param_cv.test_set);
        
        for j = 1:Y_C_test.getNumberOfExamples()
            % Get the inchis of all the candidate in the set for test
            % example j
            inchis_c = Y_C_test.getCandidateSet (j, 0, 'id');
            if (~ isnan (inchis_c)) ; continue ; end % if
            
            [~ , IX] = sort (scores{j}, 'descend');
            ranks(j) = find (strcmp (inchis_c(IX), inchis_test{j}));
        end % for
    end % for
    
    assert (sum (isnan (ranks)) == sum (isnan (mf_corres)), ...
        'The number of the examples without candidate set and number of test examples without rank (due to the absence if a candidate set) must be equal.');
    
    maxCandNum = arrayfun (@(idx) Y_C.getCandidateSet (idx, 0, 'num'), 1:Y_C.getNumberOfExamples());
    rankPerc = getRankPerc (ranks, maxCandNum);
end % function 

function [KX_list, param] = loadInputKernelList (inputDir, param)
%% LOADINPUTKERNELLIST
    switch (lower (param.data_param.inputKernel))
        case param.data_param.availInputKernels
            disp (['Evaluation using a single kernel: ', param.data_param.inputKernel]);

            KX_list = { loadKernel(strcat (inputDir, '/input_kernels/', param.data_param.inputKernel)) };

            % If a single kernel is used, than we force the "mkl" option to
            % be 'unimkl'. The weight for the kernel will be 1.
            param.mp_iokr_param.mkl = 'unimkl';
        case {'unimkl', 'alignf'}
            disp (['Evaluation using multiple kernel learning: ', param.data_param.inputKernel]);

            KX_list = cellfun (@(kernelName) loadKernel (strcat (inputDir, '/input_kernels/', kernelName)), ...
                param.data_param.availInputKernels, 'UniformOutput', false);

            param.mp_iokr_param.mkl = lower (param.data_param.inputKernel);
        otherwise
            error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', ...
                '%s: Not a valid input kernel. See MP_IOKR_DEFAULTS for a list of available kernels.', ...
                param.data_param.inputKernel);
    end % switch
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

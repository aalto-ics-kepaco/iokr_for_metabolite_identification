function IOKR_MP_reverse_feat_approximation_error (inputDir, outputDir, param)
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
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param'});    
    
    %% Load data 
    % ... input-kernels for the training examples
    if (param.debug_param.isDebugMode)
        % Lets reduce the size of the UNIMKL a bit. This is useful if we
        % want to debug the 'separate' kernel combination.
%         param.data_param.availInputKernels = {'NSF', 'PPKR', 'CPJ', 'CEC'};
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
%         param.opt_param = struct (   ...
%             'val_gamma',   [0.5, 1], ...
%             'val_lambda',  [0.5, 1], ...
%             'nOuterFolds', 10,       ...
%             'nInnerFolds', 2);           
        
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
    
    % ... pre-defined cross-validation folds to prevent to train using
    % molecular structures which are contained in the test
    ind_fold = load (strcat (inputDir, '/cv_ind.txt'));
    if (param.debug_param.isDebugMode)
        ind_fold = ind_fold(param.debug_param.debug_set);
    end % if  
       
    cv_param = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold));
    param.data_param.cv_param = cv_param;
    param.data_param.cv       = getCVIndices (param.data_param.cv_param);

    %% Evaluate the performance using 10-fold cv
    approx_error   = NaN (size (Y, 2), 1);
    approx_error_n = NaN (size (Y, 2), 1);
    mkl_weights  = NaN (numel (KX_list), param.data_param.cv.outer.NumTestSets);
    gamma_opt    = NaN (numel (KX_list), param.data_param.cv.outer.NumTestSets);  

    for foldIdx = 1:param.data_param.cv.outer.NumTestSets
        fprintf ('Outer fold: %d/%d\n', foldIdx, param.data_param.cv.outer.NumTestSets);
              
        data_param_fold = struct();
        data_param_fold.train_set = training_my (param.data_param.cv.outer, foldIdx);
        data_param_fold.test_set  = test_my (param.data_param.cv.outer, foldIdx);
        
        [approx_error(data_param_fold.test_set), mkl_weights(:, foldIdx), gamma_opt(:, foldIdx)] = ...
            approximation_error_reverse_feat (KX_list, Y, ...
                param.opt_param, param.mp_iokr_param, data_param_fold, param.debug_param);
        approx_error_n(data_param_fold.test_set) = approx_error(data_param_fold.test_set) / sum (mkl_weights(:, foldIdx));
        
        break;
            
        clear data_param_fold;
    end % for
     
    %% Store results
    result = struct ('approx_error', approx_error, 'approx_error_n', approx_error_n, ...
        'mkl_weights', mkl_weights, 'gamma_opt', gamma_opt, 'opt_param', param.opt_param, ...
        'mp_iokr_param', param.mp_iokr_param); %#ok<NASGU>
    
    settingHash = DataHash (struct (                            ...
        'cv_param',        param.data_param.cv_param,           ...
        'repetition',      param.data_param.repetition,         ...
        'input_kernel',    upper(param.data_param.inputKernel), ...
        'center',          param.mp_iokr_param.center,          ...
        'rev_iokr',        param.mp_iokr_param.rev_iokr));
    save (strcat (outputDir, '/', settingHash, '.mat'), 'result', '-v7.3');
    
    disp (settingHash);
    disp (mean (approx_error));
    disp (mkl_weights);
    disp (gamma_opt);
    disp ('! Ready !');
end % function 

function [approx_error, mkl_weights, gamma_opt] = approximation_error_reverse_feat (KX_list, Y, ...
    opt_param, mp_iokr_param, data_param, debug_param)
%% APPROXIMATION_ERROR_REVERSE_FEAT Caculates the approximation error using reverse IOKR
%   The approximation error is calculated for a test set defined as a
%   subset of the provided input kernel and output features. 
%
%   INPUTS: 
%   (data)
%       KX_list         A cell-array storing the input-kernels. 
%       Y               Matrix storing the output feature-vectors column-wise.
%   (parameter)
%       opt_param       Struct containing the parameter used for the
%                       optimization.
%       * opt_param.val_gamma       Defines the grid the gamma parameter is
%                                   searched from.
%       mp_iokr_param   Struct containing the parameter used for the
%                       reverse-iokr configuration.
%       * mp_iokr_param.rev_iokr    String either "joint" or "separate"
%                                   defining, wethere for all input kernels
%                                   the same or different gammas are used.
%       * mp_iokr_param.mkl         String either "unimkl" or "alignf"
%                                   defining the applied MKL approach. 
%       * mp_iokr_param.center      Binary indicating whether the input and
%                                   output should be centered.
%       data_param      Struct containing the parameters related to the
%                       data handling e.g.: definition of test and training
%                       set. 
%       * data_param.train_set      Binary vector indicating the examples
%                                   used for training.
%       * data_param.test_set       Binary vector indicating the examples
%                                   used for testing. The approximation
%                                   error is calculated for those examples.
%       debug_param     Struct containing the parameters related to
%                       debugging perposes.
%       * debug_param.verbose       Binary indicating whether debugging
%                                   information should be printed out.
%
%   OUTPUTS:
%       error           Vector with the length equal to the number of test
%                       examples. For each test example the approximation
%                       error is provided.
%       mkl_weights     Vector with the length equal to the number of
%                       input-kernels. For each input-kernel the MKL weight
%                       is provided, which has been calculated using the
%                       methode specified in MP_IOKR_PARAM.MKL.
%       gamma_opt       Depeding on MP_IOKR_PARAM.REV_IOKR either a single
%                       scalar value or a vector with the lenght equal to
%                       the number of input-kernels. The value(s)
%                       correspond to the optimal gamma(s) selected for the
%                       input-feature approximation.
    %% Create some stopwatches
    if (debug_param.verbose)
        sw_mkl_weights    = StopWatch ('MKL-weights');
        sw_center_norm    = StopWatch ('Center, normalize and MKL');
        sw_select_param   = StopWatch ('Select param reverse IOKR');
        sw_modify_kx_list = StopWatch ('Modify KX-list');
        sw_train_rev_iokr = StopWatch ('Train reverse IOKR');
        sw_calc_error     = StopWatch ('Calculate error');
    end % if

    %% Define training and test sets
    train_set = find (data_param.train_set);
    test_set  = find (data_param.test_set);
    
    n_test = numel (test_set);

    Y_train = Y(:, train_set);

    %% Learning kernel weights with Multiple Kernel Learning  
    KX_list_train = cellfun(@(x) x(train_set, train_set), KX_list, 'UniformOutput', false);
    
    if (debug_param.verbose)
        sw_mkl_weights.start();
    end % if
    
    mkl_weights = mkl_weight (mp_iokr_param.mkl, KX_list_train, normmat (Y_train' * Y_train));
    
    if (debug_param.verbose)
        sw_mkl_weights.stop();
        sw_mkl_weights.showAvgTime();
    end % if
    
    clear KX_list_train; 
    
    %% Centering, normalization and MKL of the input-kernels
    if (debug_param.verbose)
        sw_center_norm.start();
    end % if
    
    switch mp_iokr_param.rev_iokr
        case 'joint' 
            n_kx = 1;
            % Computation of the combined input kernel
            [KX_train_combined, KX_train_test_combined, KX_test_combined] = ...
                mkl_combine_train_test(KX_list, train_set, test_set, mkl_weights, mp_iokr_param.center);
            KX_train_list = {KX_train_combined};
            KX_train_test_list = {KX_train_test_combined};
            KX_test_list = {KX_test_combined};
            clear KX_train_combined KX_train_test_combined KX_test_combined;
        case 'separate'
            n_kx = length(KX_list);
            KX_train_list = cell(n_kx,1);
            for k = 1:n_kx
                [KX_train_list{k}, ~] = input_kernel_center_norm(KX_list{k}, train_set, test_set, mp_iokr_param.center);
                KX_train_list{k} = mkl_weights(k) * KX_train_list{k};
            end % for
        otherwise
            error ('MP_IOKR_reverse_feat:approximation_error_reverse_feat:InvalidInput', ...
                '%s is not a valid value for MP_IOKR_PARAM.REV_IOKR', mp_iokr_param.rev_iokr);
    end % switch
    
    if (debug_param.verbose)
        sw_center_norm.stop();
        sw_center_norm.showAvgTime();
    end % if
    
    %% Training output feature vectors
    mean_Y_train = mean (Y_train, 2);
    Psi_train = norma (Y_train, mean_Y_train, mp_iokr_param.center);    
    
    %% Selection of the regularization parameter(s) of reverse IOKR 
    if (debug_param.verbose)
        sw_select_param.start();
    end % if
    
    gamma_opt = Select_param_reverse_IOKR(KX_train_list, Psi_train, opt_param.val_gamma);
    
    if (debug_param.verbose)
        sw_select_param.stop();
        sw_select_param.showAvgTime();
    end % if
    
    % Modify the KX_list in order to combine the kernel belonging to unique
    % gamma values.
    
    if (debug_param.verbose)
        sw_modify_kx_list.start();
    end % if
    
    if (strcmp (mp_iokr_param.rev_iokr, 'separate'))
        gamma_opt_u = unique (gamma_opt);
        n_kx = numel (gamma_opt_u);
        
        KX_train_list      = cell (n_kx, 1);
        KX_train_test_list = cell (n_kx, 1);
        KX_test_list       = cell (n_kx, 1);
        
        k = 1;
        for val_gamma_u = gamma_opt_u'
            k_selec = (gamma_opt == val_gamma_u);
            [KX_train_list{k}, KX_train_test_list{k}, KX_test_list{k}] = mkl_combine_train_test (KX_list(k_selec), ...
                train_set, test_set, mkl_weights(k_selec), mp_iokr_param.center);
            k = k + 1;
        end % if
        fprintf ('Length of new kernel-list: %d\n', n_kx);
    end % if   
    
    if (debug_param.verbose)
        sw_modify_kx_list.stop();
        sw_modify_kx_list.showAvgTime();
    end % if    
    
    %% Training the reverse IOKR model
    if (debug_param.verbose)
        sw_train_rev_iokr.start();
    end % if
    
    M = Train_reverse_IOKR_feat(Psi_train, unique (gamma_opt));
    figure; imagesc (M);
    
    if (debug_param.verbose)
        sw_train_rev_iokr.stop();
        sw_train_rev_iokr.showAvgTime();
    end % if 
    %% Calculate error
    if (debug_param.verbose)
        sw_calc_error.start();
    end % if
    
    Psi_test = norma (Y(:, test_set), mean_Y_train, mp_iokr_param.center);
    
    I_K = repmat (eye (n_test), [n_kx, 1]);
    
%     mse = trace (cell2mat (KX_test_list') * I_K)                              ...
%         - 2 * trace (Psi_test' * M' * blkdiag (KX_train_test_list{:}) * I_K) ...
%         + trace (Psi_test' * M' * blkdiag (KX_train_list{:}) * M * Psi_test);
%     mse = mse / n_test; 
    
    approx_error = diag (cell2mat (KX_test_list') * I_K)                    ...
        - 2 * diag (Psi_test' * M' * blkdiag (KX_train_test_list{:}) * I_K) ...
        + diag (Psi_test' * M' * blkdiag (KX_train_list{:}) * M * Psi_test);
    
%     assert ((mean (approx_error) - mse) < 1e-14);
    
    if (debug_param.verbose)
        sw_calc_error.stop();
        sw_calc_error.showAvgTime();
    end % if 
end % function 
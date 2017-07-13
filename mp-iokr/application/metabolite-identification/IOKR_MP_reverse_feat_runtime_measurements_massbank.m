function IOKR_MP_reverse_feat_runtime_measurements_massbank (inputDir, outputDir, param)
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
    swc_input_kernel_processing = StopWatchCPUTime ('Input kernel processing (CPU)');
    swc_load_stats              = StopWatchCPUTime ('Load statistics (CPU)');
    swc_train_reverse           = StopWatchCPUTime ('Train model reverse (CPU)');
    swc_train                   = StopWatchCPUTime ('Train model (CPU)');
    swc_pre_image               = StopWatchCPUTime ('Pre-image (CPU)');
    
    sw_input_kernel_processing  = StopWatch ('Input kernel processing (WALL)');
    sw_load_stats               = StopWatch ('Load statistics (WALL)');
    sw_train_reverse            = StopWatch ('Train model reverse (WALL)');
    sw_train                    = StopWatch ('Train model (WALL)');
    sw_pre_image                = StopWatch ('Pre-image (WALL)');
    
    % ... input-kernels for the training examples
    % For the time measurements we use only 5 Kernels, as we do not apply
    % the "separate-trick". 
    param.data_param.availInputKernels = {'PPKR', 'NSF'};
    %     param.data_param.availInputKernels = {'PPKR', 'NSF', 'CEC', 'CPJ', 'CPJB'};
    
    swc_input_kernel_processing.start();
    sw_input_kernel_processing.start();
    
    [KX_list, param] = loadInputKernelsIntoList (strcat (inputDir, '/input_kernels/'), param);
    if (isempty (KX_list))
        error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', ...
            'No kernel loaded.');
    end % if
    
    swc_input_kernel_processing.stop();
    sw_input_kernel_processing.stop();
    
    % Set seed in order to be able to reproduce the settings.
    % By default this is set to 'shuffle'. This means, that the current
    % time is used as seed. However, if we run the script on triton the
    % jobs might start at the same time. We therefore need to be able to
    % set the seed manualy. For example by using the slurm-job-id. This can
    % be done by the calling sbatch-file.
    rng (param.debug_param.randomSeed);
    
    % For the time-measurements we can fix the gamma and lambda parameter.
    % Furthermore we are only using one outer fold for the training. This
    % means all the data are going to be training data. We test on an
    % independent subset of 625 Massbank data.
    param.opt_param = struct ( ...
        'val_gamma',   1,      ...
        'val_lambda',  1,      ...
        'nOuterFolds', 1,      ...
        'nInnerFolds', 10);        
    
    if (param.debug_param.isDebugMode)
        fprintf ('Show me a random-number: %f\n', rand (1));
        
        % Modify the output-dir & creat it if needed
        outputDir = strcat (outputDir, '/debug/');
        if (~ exist (outputDir, 'dir'))
            if (~ mkdir (outputDir))
                error ('IOKR_MP_reverse_feat_evaluation:RuntimeError', ...
                    'Could not create debug-directory: %s', outputDir);
            end % if
        end % if
               
        n = size (KX_list{1}, 1);
        param.debug_param.debug_set = false (n, 1);        
        param.debug_param.debug_set(randsample (n, param.debug_param.n_debug_set)) = true;
        KX_list = cellfun(@(x) x(param.debug_param.debug_set, param.debug_param.debug_set), KX_list, 'UniformOutput', false);
    end % if
    n_kx = numel (KX_list);
    
    % ... fingerprints for the training examples
    Y = load (strcat (inputDir, '/fingerprints/fp.mat'));
    Y = full (Y.Y);
    if (param.debug_param.isDebugMode)
        Y = Y(:, param.debug_param.debug_set);
    end % if
    n_train = size (Y, 2);
    
    % ... identifier (inchis) for the training examples
    inchis = readtext (strcat (inputDir, '/inchi.txt'));
    if (param.debug_param.isDebugMode)
        inchis = inchis(param.debug_param.debug_set);
    end % if
    
    %% Load / Store pre-calculated statistics
    cv_param = struct ('nObservations', n_train,                        ...
                       'outer', struct ('type', 'random', 'nFolds', 1), ...
                       'inner', struct ('nFolds', param.opt_param.nInnerFolds));                        
    param.data_param.cv_param = cv_param;
    param.data_param.selection_param.inclExpCand = false;

    % NOTE: The selec_ property of Y_C will be modified according to
    %       the selection defined by PARAM.DATA_PARAM.SELECTION_PARAM.
    tic;
    matObj = getPreCalcCandStat_feat (Y, [], inchis, param, strcat (inputDir, '/pre_calculated_stats/'), true);
    fprintf ('Loading / pre-calculating of the candidate statistics took %.3fs\n', toc);

    param.data_param.cv         = matObj.cv;
    param.mp_iokr_param.center  = matObj.center;
    param.data_param.repetition = matObj.repetition;
    % matObj also contains the statistics
    param.data_param.matObj     = matObj;
    
    %% Run the training
    % We use all the data for training - 4138 GNPS data. The test set can
    % be empty.
    
    % Centering, normalization and MKL of the input-kernels
    w = repmat (1 / n_kx, [n_kx, 1]);
    
    swc_input_kernel_processing.start();
    sw_input_kernel_processing.start();
    
    switch param.mp_iokr_param.rev_iokr
        case 'joint' 
            KX_train_list = cell (1, 1);
            % Computation of the combined input kernel
            KX_train_list{1} = zeros (size (KX_list{1}));
            for k = 1:n_kx
                KX_train_list{1} = KX_train_list{1} + w(k) * ...
                    normmat ( ...
                    center (  ...
                        KX_list{k}, mean(KX_list{k}, 1), param.mp_iokr_param.center));                   
            end % for
        case 'separate'
            KX_train_list = cell(n_kx,1);
            for k = 1:n_kx
                KX_train_list{k} = normmat ( ...
                                   center (  ...
                                        KX_list{k}, mean(KX_list{k}, 1), param.mp_iokr_param.center));
                KX_train_list{k} = KX_train_list{k} * w(k);
            end
        otherwise
            error ('MP_IOKR_reverse_feat:InvalidInput', ...
                '%s is not a valid value for MP_IOKR_PARAM.REV_IOKR', param.mp_iokr_param.rev_iokr);
    end % switch 
    
    swc_input_kernel_processing.stop();
    sw_input_kernel_processing.stop();
    
    % Training output feature vectors
    Y_train = Y;
    mean_Y_train = mean (Y_train,2);
    Psi_train = norma (Y_train, mean_Y_train, param.mp_iokr_param.center);
       
    % Training the reverse IOKR model
    swc_train_reverse.start();
    sw_train_reverse.start();
    
    M = Train_reverse_IOKR_feat (Psi_train, repmat (param.opt_param.val_gamma, [numel(KX_train_list), 1]));
    
    swc_train_reverse.stop();
    swc_train_reverse.showAvgTime();
    sw_train_reverse.stop();
    sw_train_reverse.showAvgTime();

    % Training the MP-IOKR model
    swc_load_stats.start();
    sw_load_stats.start();
    
    stats = param.data_param.matObj.stats(1, 1);
    Mean_Psi_C_train = stats.Mean_Psi_C_train;
    Cov_Psi_C_train  = stats.Cov_Psi_C_train;
    
    swc_load_stats.stop();
    swc_load_stats.showAvgTime();
    sw_load_stats.stop();
    sw_load_stats.showAvgTime();
    
    swc_train.start();
    sw_train.start();
    
    C = Train_MP_IOKR_reverse_feat (KX_train_list, Psi_train, M, ...
        Mean_Psi_C_train, Cov_Psi_C_train, ...
        param.opt_param.val_lambda);
    
    swc_train.stop();
    swc_train.showAvgTime();
    sw_train.stop();
    sw_train.showAvgTime();

    %% Run the test
    swc_input_kernel_processing.start();
    sw_input_kernel_processing.start();
    
    [KX_train_test_list, param] = loadInputKernelsIntoList ( ...
        strcat (inputDir, '/input_kernels/massbank/gnps_massbank/'), param);
    [KX_test_list, param] = loadInputKernelsIntoList ( ...
        strcat (inputDir, '/input_kernels/massbank/massbank/'), param);
    
    eval_set = load (strcat (inputDir, '/input_kernels/massbank/gnps_massbank/ind_eval.txt'));
    KX_train_test_list = cellfun(@(x) x(:, eval_set), KX_train_test_list, 'UniformOutput', false);
    
    switch param.mp_iokr_param.rev_iokr
        case 'joint' 
            % Computation of the combined input kernel
            KX_train_test_combined = zeros (size (KX_train_test_list{1}));
            for k = 1:n_kx
                mean_KX_train      = mean (KX_list{k}, 1);
                mean_KX_train_test = mean (KX_train_test_list{k}, 1);
                
                KX_train_c    = center (KX_list{k}, mean(KX_list{k}, 1), param.mp_iokr_param.center);
                KX_test_c     = center (KX_test_list{k}, mean_KX_train, param.mp_iokr_param.center, ...
                    mean_KX_train_test', mean_KX_train_test);
                KX_train_test_cn = center (KX_train_test_list{k}, mean_KX_train, param.mp_iokr_param.center, ...
                    mean_KX_train', mean_KX_train_test);
                KX_train_test_cn = KX_train_test_cn ./ sqrt (diag (KX_train_c) * diag (KX_test_c)');
                
                KX_train_test_combined = KX_train_test_combined + w(k) * KX_train_test_cn;               
            end % for
            KX_train_test_list = {KX_train_test_combined};
            
            clear KX_train_test_combined KX_test_list KX_train_c KX_test_c KX_train_test_cn;
        case 'separate'
            for k = 1:n_kx
                mean_KX_train      = mean (KX_list{k}, 1);
                mean_KX_train_test = mean (KX_train_test_list{k}, 1);
                
                KX_train_c    = center (KX_list{k}, mean(KX_list{k}, 1), param.mp_iokr_param.center);
                KX_test_c     = center (KX_test_list{k}, mean_KX_train, param.mp_iokr_param.center, ...
                    mean_KX_train_test', mean_KX_train_test);
                KX_train_test_cn = center (KX_train_test_list{k}, mean_KX_train, param.mp_iokr_param.center, ...
                    mean_KX_train', mean_KX_train_test);
                KX_train_test_cn = KX_train_test_cn ./ sqrt (diag (KX_train_c) * diag (KX_test_c)');
                
                KX_train_test_list{k} = w(k) * KX_train_test_cn;
            end
            
            clear KX_test_list KX_train_c KX_test_c KX_train_test_cn;
        otherwise
            error ('MP_IOKR_reverse_feat:InvalidInput', ...
                '%s is not a valid value for MP_IOKR_PARAM.REV_IOKR', param.mp_iokr_param.rev_iokr);
    end % switch 
    
    swc_input_kernel_processing.stop();
    sw_input_kernel_processing.stop();
    swc_input_kernel_processing.showAvgTime();
    sw_input_kernel_processing.showAvgTime();
    
    % Pre-image
    swc_pre_image.start();
    sw_pre_image.start();
    
    mf_test = readtext(strcat(inputDir, '/candidates/massbank/MF.txt'));
    n_test = length(mf_test);
    inchi_pred = cell(n_test,1);
    
    Psi_pred = Prediction_MP_IOKR_reverse_feat(KX_train_test_list, C);
    
    for j = 1:n_test
        % Candidate inchi
        inchi_c = readtext(strcat(inputDir, '/candidates/massbank/inchi/inchi_',mf_test{j},'.txt'));

        % Candidate fingerprints
        load(strcat(inputDir, '/candidates/massbank/fp/fp_',mf_test{j},'.mat'));
        Y_Cj = full(fp');
        Psi_Cj = norma(Y_Cj, mean_Y_train, param.mp_iokr_param.center);

        % Computation of the score function between the test example
        % and the corresponding candidates
        score = Psi_pred(:, j)' * Psi_Cj;
        [~,IX] = max(score);
        inchi_pred{j} = inchi_c{IX};
    end
    
    swc_pre_image.stop();
    sw_pre_image.stop();
    swc_pre_image.showAvgTime();
    sw_pre_image.showAvgTime();
    
    disp ('! Ready !');
end % function 

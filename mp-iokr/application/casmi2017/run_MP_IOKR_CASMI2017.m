function run_MP_IOKR_CASMI2017 (input_dir_training, input_dir_test, output_dir, mpiokr_param, challenge_param)
%======================================================
% DESCRIPTION:
% Script for running MP-IOKR on a small test-dataset containing ~260
% MS2-spectra. 
%
% INPUTS:
% inputDir:     directory in which the data are contained
% result_dir:   directory in which the results will be saved
% cand:         candidate set
%
% NOTES:  
%   - Evaluate only on those examples with a candidate-set size smaller
%     than X:
%   eval = find (arrayfun (@(x) cand(x).num < 3000, mf_corres));
%
%======================================================

    %--------------------------------------------------------------
    % Set up parameters
    %--------------------------------------------------------------
    %param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
    %    {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param', 'ky_param'});
      
%     param.opt_param.nInnerFolds = 2;
    
    %--------------------------------------------------------------
    % Load and prepare data
    %--------------------------------------------------------------

    % inchi keys, molecular formulas, fingerprints
    load ([input_dir_training, '/compounds/', challenge_param.ion_mode, '/compound_info.mat'], 'dt_inchi_mf_fp');
    inchi = dt_inchi_mf_fp.inchi_key_1; 
    
    % Extract fingerprints
    switch (challenge_param.fp_set)
        case 'full'
            Y_train = full (dt_inchi_mf_fp.fp_full)';
        case 'masked'
            Y_train = full (dt_inchi_mf_fp.fp_masked)';
    end % switch
    
    % Load fingerprint masks
    switch (challenge_param.fp_set)
        case 'full'
            fp_mask = 'ALL';
        case 'masked'
            fp_mask_fn = [input_dir_training, '/fingerprints/', challenge_param.ion_mod, '/fingerprints.mask'];
            fp_mask = load_fingerprint_mask (fp_mask_fn);
    end % switch
    
    % Candidates description
    tic;
    load ([input_dir_training, '/candidates/', challenge_param.ion_mode, '/cand_prepared.mat'], 'cand');
    toc;
    mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, cand);
    Y_C_train = CandidateSets (DataHandle (cand), mf_corres, 'ALL', fp_mask);  
    assert (Y_C_train.getNumberOfExamples() == size (Y_train, 2))
    
    % Candidate selection
    selec = getCandidateSelection (Y_C_train, inchi, mpiokr_param.data_param.selection_param);      
    Y_C_train.setSelectionsOfCandidateSets (selec);
    
    % Input kernels
    kernel_dir = [input_dir_training, '/kernels/', challenge_param.ion_mode, '/'];
    kernel_files = dir ([kernel_dir, '/*.mat']);
    mpiokr_param.data_param.availInputKernels = arrayfun (@(file) basename (file.name), ...
        kernel_files, 'UniformOutput', false);
    mpiokr_param.data_param.inputKernel = '__ALL__';
    KX_list_train = loadInputKernelsIntoList (kernel_dir, mpiokr_param, '.mat');
    
    %--------------------------------------------------------------
    % Train the model using all the training data
    %--------------------------------------------------------------
    % Initiate the output kernel approximation in case
    if (strcmp (mpiokr_param.ky_param.type, 'gaussian')) && (strcmp (mpiokr_param.ky_param.representation, 'feature'))
        % RANDOM FOURIER FEATURES
        assert (ismember ('rff_dimension', fieldnames (mpiokr_param.ky_param)), ...
            'The random fourier feature dimension must be specified.');
        
        mpiokr_param.ky_param.rff = RandomFourierFeatures (size (Y_train, 1), mpiokr_param.ky_param.rff_dimension);
    end % if
    
    train_model = Train_MPIOKR (KX_list_train, Y_train, Y_C_train, ...
            mpiokr_param.ky_param, mpiokr_param.mp_iokr_param, mpiokr_param.opt_param, mpiokr_param.debug_param);
    
    filename = [output_dir, '/', challenge_param.ion_mode, '/train_model_mpiokr_mkl=', iokr_param.iokr_param.mkl, ...
       '_kernel=', mpiokr_param.ky_param.type, '_representation=', mpiokr_param.ky_param.representation, ...
       '_base=', mpiokr_param.ky_param.base_kernel, '_', mpiokr_param.ky_param.param_selection, ...
       '_strategy=', mpiokr_param.data_param.selection_param.strategy, ...
       '_ion_mode=', challenge_param.ion_mode, ...
       '_fp_set=', challenge_param.fp_set, '.mat'];
    save (filename, 'train_model', '-v7.3');
        
    %--------------------------------------------------------------
    % Get the scoring for each challenge
    %--------------------------------------------------------------
    switch challenge_param.ion_mode
        case 'positive'
            challenge_ids = 132:243;
        case 'negative'
            challenge_ids = 46:131;
    end % switch

    for challenge_idx = challenge_ids
        fprintf ('Calculate scores for challenge %d\n', challenge_idx);
        
        % Load molecular candidates
        challenge_cand_fn = sprintf ( ...
            '%s/%s/candidates/%s/candidates-challenge-%03d-inchi-level1-nonredundant.fpt.mat', ...
            input_dir_test, challenge_param.ion_mode, challenge_param.fp_set, challenge_idx);
        fprintf ('Filename candidates: %s\n', challenge_cand_fn);
        
        tmp         = load (challenge_cand_fn);
        tmp.cand.id = tmp.cand.id';
        Y_C_test    = CandidateSets (DataHandle (tmp.cand), 1);
        
        inchis_test = Y_C_test.getCandidateSet (1, false, 'id');
        
        result_test = table (inchis_test', 'VariableNames', {'inchi'});
        
        for mf_cand_id = 1:1
            % Load the train vs. test and test vs. test kernels
            challenge_KX_train_test_fn = dir (sprintf ( ...
                '%s/kernels_training_test/challenge-%03d-msms.mgf_01_*.kernels', ...
                input_dir_test, challenge_idx));
            challenge_KX_test_fn       = dir (sprintf ( ...
                '%s/kernels_test_test/challenge-%03d-msms.mgf_01_*.test_kernels', ...
                input_dir_test, challenge_idx));
            
            if (isempty (challenge_KX_train_test_fn))
                break;
            end % if

            challenge_KX_train_test_fn = challenge_KX_train_test_fn.name;
            challenge_KX_test_fn       = challenge_KX_test_fn.name;

            fprintf ('Filename KX_train_test: %s\n', challenge_KX_train_test_fn)
            fprintf ('Filename KX_test: %s\n',       challenge_KX_test_fn)         
            
            [KX_list_train_test, KX_names] = read_challenge_kernel ( ...
            [input_dir_test, '/kernels_training_test/', challenge_KX_train_test_fn]);
            [~, locb]                      = ismember (upper (KX_names), mpiokr_param.data_param.availInputKernels);
            KX_list_train_test             = KX_list_train_test(locb);
        
            [KX_list_test, KX_names] = read_challenge_kernel ( ...
            [input_dir_test, '/kernels_test_test/', challenge_KX_test_fn]);
            [~, locb]                = ismember (upper (KX_names), mpiokr_param.data_param.availInputKernels);
            KX_list_test             = KX_list_test(locb);
            
            % Calcualte and write out scores
            scores_test = Test_MPIOKR (KX_list_train_test, KX_list_test, train_model, ...
                Y_train, Y_C_train, Y_C_test, mpiokr_param.mp_iokr_param, mpiokr_param.mp_iokr_param.center, ...
                mpiokr_param.debug_param);
            
            result_test = [result_test, ...
                table(scores_test{1}', 'VariableNames', {sprintf('score_%02d', mf_cand_id)})];
        end % for
        
        writetable (result_test, sprintf ('%s/%s/mpiokr-category4-%03d.txt', ...
            output_dir, challenge_param.ion_mode, challenge_idx), ...
            'Delimiter', '\t', 'WriteVariableNames', true);
    end % for
end
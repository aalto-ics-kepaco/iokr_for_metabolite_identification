function run_IOKR_CASMI2017 (input_dir_training, input_dir_test, output_dir, ...
    iokr_param, challenge_param)
    %% RUN_IOKR_CASMI2017 Scoring script for the CASMI2017 challenge spectra
    %   INPUTS:
    %   input_dir_training: Base-directory containing the training data.
    %   input_dir_test: Base-directory containing the challenge data.
    %   output_dir: Base-directory to store the trained models and scores.
    %   iokr_param: Structure containing the parameters for the IOKR algorithm.
    %          See MP_IOKR_Defaults for details.
    %   challenge_param: Structure containing the parameters regarding the 
    %                    challenge data:
    %                    - ion_mode: String {'positive', 'negative'}
    %                    - fp_set: String {'full', 'masked'}
    
    %--------------------------------------------------------------
    % Load and prepare data
    %--------------------------------------------------------------
    % Fingerprints
    load ([input_dir_training, '/compounds/', challenge_param.ion_mode, '/compound_info.mat'], 'dt_inchi_mf_fp');
    switch (challenge_param.fp_set)
        case 'full'
            Y_train = full (dt_inchi_mf_fp.fp_full)';
        case 'masked'
            Y_train = full (dt_inchi_mf_fp.fp_masked)';
    end % switch
            
    % Input kernels
    kernel_dir = [input_dir_training, '/kernels/', challenge_param.ion_mode, '/'];
    kernel_files = dir ([kernel_dir, '/*.mat']);
    iokr_param.data_param.availInputKernels = arrayfun (@(file) basename (file.name), ...
        kernel_files, 'UniformOutput', false);
    iokr_param.data_param.inputKernel = '__ALL__';
    KX_list_train = loadInputKernelsIntoList (kernel_dir, iokr_param, '.mat');
    
    %--------------------------------------------------------------
    % Train the model using all the training data
    %--------------------------------------------------------------
    train_model = Train_IOKR (KX_list_train, Y_train, ...
            iokr_param.ky_param, iokr_param.opt_param, iokr_param.iokr_param);
    
    filename = [output_dir, '/', challenge_param.ion_mode, '/train_model_iokr_mkl=', iokr_param.iokr_param.mkl, ...
       '_kernel=', iokr_param.ky_param.type, ...
       '_base=', iokr_param.ky_param.base_kernel, '_', iokr_param.ky_param.param_selection, ...
       '_ion_mode=', challenge_param.ion_mode, ...
       '_fp_set=', challenge_param.fp_set, '.mat'];
    save (filename, 'train_model', '-v7.3');
%     load (filename);    
     
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
        
        for mf_cand_id = 1:5
            % Load the train vs. test and test vs. test kernels
            challenge_KX_train_test_fn = dir (sprintf ( ...
                '%s/%s/kernels_training_test/challenge-%03d-msms.mgf_%02d_*.kernels', ...
                input_dir_test, challenge_param.ion_mode, challenge_idx, mf_cand_id));
            challenge_KX_test_fn       = dir (sprintf ( ...
                '%s/%s/kernels_test_test/challenge-%03d-msms.mgf_%02d_*.test_kernels', ...
                input_dir_test, challenge_param.ion_mode, challenge_idx, mf_cand_id));

            if (isempty (challenge_KX_train_test_fn))
                break;
            end % if
            
            challenge_KX_train_test_fn = challenge_KX_train_test_fn.name;
            challenge_KX_test_fn       = challenge_KX_test_fn.name;

            fprintf ('Filename KX_train_test: %s\n', challenge_KX_train_test_fn)
            fprintf ('Filename KX_test: %s\n',       challenge_KX_test_fn)

            [KX_list_train_test, KX_names] = read_challenge_kernel ( ...
                [input_dir_test, '/', challenge_param.ion_mode, '/kernels_training_test/', challenge_KX_train_test_fn]);
            [~, locb]                      = ismember (upper (KX_names), iokr_param.data_param.availInputKernels);
            KX_list_train_test             = KX_list_train_test(locb);

            [KX_list_test, KX_names] = read_challenge_kernel ( ...
                [input_dir_test, '/', challenge_param.ion_mode, '/kernels_test_test/', challenge_KX_test_fn]);
            [~, locb]                = ismember (upper (KX_names), iokr_param.data_param.availInputKernels);
            KX_list_test             = KX_list_test(locb);

            % Calcualte and write out scores
            scores_test = Test_IOKR (KX_list_train_test, KX_list_test, train_model, ...
                Y_train, Y_C_test, iokr_param.iokr_param.center);

            result_test = [result_test, table(scores_test{1}'+abs(min (scores_test{1})), ... 
                'VariableNames', {sprintf('score_%02d', mf_cand_id)})];
        end % for
            
        writetable (result_test, sprintf ('%s/%s/iokr-category4-%03d.txt', ...
            output_dir, challenge_param.ion_mode, challenge_idx), ...
            'Delimiter', '\t', 'WriteVariableNames', false);
    end % for
end
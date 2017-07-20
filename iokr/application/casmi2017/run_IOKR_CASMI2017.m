function run_IOKR_CASMI2017 (input_dir_training, input_dir_test, output_dir, param)
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
   % param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
     %   {'debug_param', 'opt_param', 'iokr_param', 'data_param', 'ky_param'});
%     
   % param.iokr_param.model_representation = 'Chol_decomp_of_C';
    %--------------------------------------------------------------
    % Load and prepare data
    %--------------------------------------------------------------

    % inchi keys, molecular formulas, fingerprints
    load ([input_dir_training '/compounds/compound_info.mat'], 'dt_inchi_mf_fp');
    
    % Extract fingerprints
    Y_train = full (dt_inchi_mf_fp.fp_masked)';
    %param.ky_param.representation  = 'feature';
    %param.ky_param.type            = 'linear';
    %param.ky_param.base_kernel     = 'linear';
    %param.ky_param.param_selection = 'cv';
    
    % Input kernels
    kernel_files = dir ([input_dir_training '/kernels/*.mat']);
    param.data_param.availInputKernels = arrayfun (@(file) basename (file.name), ...
        kernel_files, 'UniformOutput', false);
    %param.data_param.availInputKernels = param.data_param.availInputKernels(1:2);
    param.data_param.inputKernel = param.iokr_param.mkl;
    KX_list_train = loadInputKernelsIntoList ([input_dir_training, '/kernels/'], param, '.mat');
    
    %--------------------------------------------------------------
    % Train the model using all the training data
    %--------------------------------------------------------------
    train_model = Train_IOKR (KX_list_train, Y_train, ...
            param.ky_param, param.opt_param, param.iokr_param);
%    save ([input_dir_training, '/train_model.mat'], 'train_model', '-v7.3');
    
    filename = [input_dir_training, '/train_model_iokr_mkl=', param.iokr_param.mkl, ...
       '_kernel=', param.ky_param.type, ...
       '_base=', param.ky_param.base_kernel, '_', param.ky_param.param_selection, '.mat'];
    save (filename, 'train_model', '-v7.3');
        
    %--------------------------------------------------------------
    % Get the scoring for each challenge
    %--------------------------------------------------------------
    % TODO: FOR THE MOMENT WE ONLY CONSIDER THE MOST CANDIDATE MOLECULAR 
    %       FORMULA WITH THE HIGHEST SCORE.
    
    challenge_ids = 46:243;
    for challenge_idx = challenge_ids
        fprintf ('Calculate scores for challenge %d\n', challenge_idx);
        
        % Load the train vs. test and test vs. test kernels
        challenge_KX_train_test_fn = dir (sprintf ( ...
            '%s/kernels_training_test/challenge-%03d-msms.mgf_01_*.kernels', ...
            input_dir_test, challenge_idx));
        challenge_KX_test_fn       = dir (sprintf ( ...
            '%s/kernels_test_test/challenge-%03d-msms.mgf_01_*.test_kernels', ...
            input_dir_test, challenge_idx));
        
        challenge_KX_train_test_fn = challenge_KX_train_test_fn.name;
        challenge_KX_test_fn       = challenge_KX_test_fn.name;
        
        fprintf ('Filename KX_train_test: %s\n', challenge_KX_train_test_fn)
        fprintf ('Filename KX_test: %s\n',       challenge_KX_test_fn)
        
        [KX_list_train_test, KX_names] = read_challenge_kernel ( ...
            [input_dir_test, '/kernels_training_test/', challenge_KX_train_test_fn]);
        [~, locb]                      = ismember (upper (KX_names), param.data_param.availInputKernels);
        KX_list_train_test             = KX_list_train_test(locb);
        
        [KX_list_test, KX_names] = read_challenge_kernel ( ...
            [input_dir_test, '/kernels_test_test/', challenge_KX_test_fn]);
        [~, locb]                = ismember (upper (KX_names), param.data_param.availInputKernels);
        KX_list_test             = KX_list_test(locb);
        
%         KX_list_train_test = cellfun (@(x) x(1:260, 1), KX_list_train_test, ...
%             'UniformOutput', false);
%         
        % Load molecular candidates
        challenge_cand_fn = sprintf ( ...
            '%s/candidates-challenge-%03d-inchi-level1-nonredundant.fpt.mat', ...
            input_dir_test, challenge_idx);
        fprintf ('Filename candidates: %s\n', challenge_cand_fn);
        
        tmp         = load (challenge_cand_fn);
        tmp.cand.id = tmp.cand.id';
%         tmp.cand.data = tmp.cand.data(1:709,:);
        Y_C_test    = CandidateSets (DataHandle (tmp.cand), 1);
        
        % Calcualte and write out scores
        scores_test = Test_IOKR (KX_list_train_test, KX_list_test, train_model, ...
            Y_train, Y_C_test, param.iokr_param.center);
        
        inchis_test = Y_C_test.getCandidateSet (1, false, 'id');
        
        result_test = table (inchis_test', scores_test{1}' + abs (min (scores_test{1})), ... 
            'VariableNames', {'inchi', 'score'});
        writetable (result_test, sprintf ('%s/iokr-4-%03d.txt', output_dir, challenge_idx), ...
            'Delimiter', '\t', 'WriteVariableNames', false);
    end % for
end
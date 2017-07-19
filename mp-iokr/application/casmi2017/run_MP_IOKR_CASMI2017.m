function run_MP_IOKR_CASMI2017 (input_dir_training, input_dir_test, output_dir)
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
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
        {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param', 'ky_param'});
      
    %--------------------------------------------------------------
    % Load and prepare data
    %--------------------------------------------------------------

    % inchi keys, molecular formulas, fingerprints
    load ([input_dir_training '/compound_info.mat'], 'dt_inchi_mf_fp');
    inchi = dt_inchi_mf_fp.inchi_key_1; 
    
    % Extract fingerprints
    Y_train = full (dt_inchi_mf_fp.fp_masked)';
    [~,n] = size(Y_train);
    param.ky_param.representation  = 'feature';
    param.ky_param.type            = 'linear';
    param.ky_param.base_kernel     = 'linear';
    param.ky_param.param_selection = 'cv';
   
    % Cross-validation
    cv = getCVIndices (struct ('nObservations', n, ...
        'outer', struct ('type', 'random', 'nFolds', n_folds)));

    % Candidates description
    mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, cand);
    Y_C_train = CandidateSets (DataHandle (cand), mf_corres);  
    assert (Y_C_train.getNumberOfExamples() == size (Y_train, 2))
    
    % Candidate selection
    % TODO: WE USE ALL THE CANDIDATES FOR THE MOMENT
    
    % Input kernels
    kernel_files = dir ([input_dir_training '/kernels/*.txt']);
    param.data_param.availInputKernels = arrayfun (@(file) basename (file.name), ...
        kernel_files, 'UniformOutput', false);
    param.data_param.inputKernel = 'unimkl';
    KX_list_train = loadInputKernelsIntoList ([input_dir_training, '/kernels/'], param, '.txt');
    
    %--------------------------------------------------------------
    % Train the model using all the training data
    %--------------------------------------------------------------
    train_model = Train_MPIOKR (KX_list_train, Y_train, Y_C_train, ...
            param.ky_param, param.mp_iokr_param, param.opt_param, param.debug_param);
    save ([input_dir_training, '/train_model.mat'], 'train_model', '-v7.3');
    
    %    filename = [output_dir, '/mpiokr/', 'rank_mkl=', param.mp_iokr_param.mkl, ...
    %    '_kernel=', param.ky_param.type, ...
    %    '_base=', param.ky_param.base_kernel, ...
    %    '_', param.ky_param.param_selection, ...
    %    '_strategy=', param.data_param.selection_param.strategy, ...
    %    '_inclCandExp', num2str(param.data_param.selection_param.inclCandExp)];
    %save(filename,'rank_perc_100','-ascii');
        
    %--------------------------------------------------------------
    % Get the scoring for each challenge
    %--------------------------------------------------------------
    % TODO: FOR THE MOMENT WE ONLY CONSIDER THE MOST CANDIDATE MOLECULAR 
    %       FORMULA WITH THE HIGHEST SCORE.
    
    challenge_ids = 46:243;
    for challenge_idx = challenge_ids'
        fprintf ('Calculate scores for challenge %d\n', challenge_idx);
        
        % Load the train vs. test and test vs. test kernels
        challenge_KX_train_test_fn = dir (sprintf ( ...
            '%s/kernels_train_test/challenge-%03d-msms.mfg-_01_*.kernels', ...
            input_dir_test, challenge_idx));
        challenge_KX_test_fn       = dir (sprintf ( ...
            '%s/kernels_test_test/challenge-%03d-msms.mfg-_01_*.test_kernels', ...
            input_dir_test, challenge_idx));
        
        challenge_KX_train_test_fn = challenge_KX_train_test_fn.name;
        challenge_KX_test_fn       = challenge_KX_test_fn.name;
        
        fprintf ('Filename KX_train_test: %s\n', challenge_KX_train_test_fn)
        fprintf ('Filename KX_test: %s\n',       challenge_KX_test_fn)
        
        KX_list_train_test = read_challenge_kernel (challenge_KX_train_test_fn);
        KX_list_test       = read_challenge_kernel (challenge_KX_test_fn);
        
        % Load molecular candidates
        challenge_cand_fn = sprintf ( ...
            '%s/candidates_nonredundant/candidates-challenge-%03d-inchi-level1-nonredundant.fpt.mat', ...
            input_dir_test, challenge_idx);
        fprintf ('Filename candidates: %s\n', challenge_cand_fn);
        
        tmp      = load (challenge_cand_fn);
        Y_C_test = CandidateSets (tmp.cand, 1);
        
        % Calcualte and write out scores
        scores_test = Test_MPIOKR (KX_list_train_test, KX_list_test, train_model, ...
            Y_train, Y_C_train, Y_C_test, param.mp_iokr_param, param.mp_iokr_param.center, param.debug_param);
        
        inchis_test = Y_C_test.getCandidateSet (1, false, 'id');
        
        result_test = table (inchis_test, scores_test, 'VariableNames', {'inchi', 'score'});
        writetable (result_test, sprintf ('%s/mpiokr-4-%03d.txt', output_dir, challenge_idx), ...
            'Delimiter', '\t', 'WriteVariableNames', false);
    end % for
end
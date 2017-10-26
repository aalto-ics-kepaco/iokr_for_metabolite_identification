function train_models_for_csifingerid_cv (wdir)
%======================================================
% DESCRIPTION:
% Script for running IOKR
%
% INPUTS:
% inputDir:     directory in which the data are contained
% result_dir:   directory in which the results will be saved
%
%======================================================
    outputDir = [wdir, '/models_iokr/'];
    if (~ exist (outputDir, 'dir'))
        mkdir (outputDir);
    end % if

    %--------------------------------------------------------------
    % Load Data
    %--------------------------------------------------------------
    % inchi keys, molecular formulas, fingerprints
    if (~ exist ([wdir, '/spectra/compound_info.mat'], 'file'))
        prepare_compound_info_for_iokr (wdir)
    end % if
    
    load ([wdir, '/spectra/compound_info.mat'], 'dt_inchi_mf_fp'); 
    Y = full (dt_inchi_mf_fp.fp_full)';
    
    % Fingerprint mask 
    fp_mask = load_fingerprint_mask ([wdir, '/fingerprints/fingerprints.mask']);
    Y = Y(fp_mask, :);
    
    % Run IOKR with default parameters
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct ( ...
        'ky_param', struct ( ...    
            'representation',  'kernel',   ...
            'type',            'gaussian', ...
            'base_kernel',     'tanimoto', ...
            'param_selection', 'entropy'), ...
        'iokr_param', struct ( ...
            'mkl', 'alignf')), ...
        {'opt_param', 'iokr_param', 'data_param', 'ky_param'});       
    
    % Find available input kernels
    kernel_files = dir ([wdir '/kernels/*.txt']);
    if isempty (kernel_files)
        error ('No kernels available. Did you change their file-extensions to ".txt"?');
    end % if
    param.data_param.availInputKernels = arrayfun (@(file) basename (file.name), ...
        kernel_files, 'UniformOutput', false);
    
    % Convert kernel matrices to matfiles, if those do not exist yet.
    for input_kernel = param.data_param.availInputKernels'
        if (~ exist ([wdir, '/kernels/', input_kernel{1}, '.mat'], 'file'))
            convert_kernel_file_to_mat ([wdir, '/kernels/'], input_kernel{1});
        end % for
    end % for  
    KX_list = loadInputKernelsIntoList ([wdir '/kernels/'], param, '.mat');
    
    % Pre-defined cross-validation folds
    dt_cv_folds = readtable ([wdir, '/crossvalidation_folds.txt'], ...
        'Delimiter', '\t', 'ReadVariableNames', false);
    dt_cv_folds.Properties.VariableNames = {'mol_id', 'inchi_key_1', 'cv'};
    dt_inchi_mf_fp = join (dt_inchi_mf_fp, dt_cv_folds, 'Keys', 'mol_id');
    
    %--------------------------------------------------------------
    % Model parameter estimation using nested cross-validation
    %--------------------------------------------------------------
    disp ('Model training')
    
    n_folds = max (dt_inchi_mf_fp.cv);
    for cv_idx = 1:n_folds
        filename = [outputDir, 'model_', MP_IOKR_Defaults.param2str(param), ...
            '_cv=', num2str(cv_idx), '.mat'];
        
        if (~ exist (filename, 'file'))
            train_set = find (dt_inchi_mf_fp.cv ~= cv_idx);
            
            KX_list_train = cellfun(@(x) x(train_set, train_set), KX_list, ...
                'UniformOutput', false);
            Y_train = Y(:, train_set);
            
            train_model = Train_IOKR (KX_list_train, Y_train, ...
                param.ky_param, param.opt_param, param.iokr_param);

            switch train_model.KY_par.representation
                case 'feature'
                    [~, process_output] = output_feature_preprocessing_train (Y_train, ...
                        train_model.ker_center);
                case 'kernel'
                    [~, process_output] = output_kernel_preprocessing_train (Y_train, ...
                        train_model.KY_par, train_model.ker_center);
            end % switch 
            train_model.process_output = process_output;
            train_model.kernel_names = param.data_param.availInputKernels;

            save (filename, 'train_model', '-v7.3');
            
    %--------------------------------------------------------------
    % Save the model: As binary for csi-fingerid 
    %--------------------------------------------------------------
            outputDir_binary = [outputDir, '/binary/cv=', num2str(cv_idx)];
            if (~ exist (outputDir_binary, 'dir'))
                mkdir (outputDir_binary);
            end % if
            write_out_iokr_model_as_binary_files (train_model, outputDir_binary);
        end % if
    end % for
    
    %--------------------------------------------------------------
    % Prediction & scoring 
    %--------------------------------------------------------------
    disp ('Prediction & Scoring')
    error ('Not implemented yet.')
    
    mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, cand);
    load ([wdir, '/cand.mat'], 'cand');
    Y_C = CandidateSets (DataHandle (cand), mf_corres, 'ALL', fp_mask);
    
    scores = cell (Y_C.getNumberOfExamples(), 1);
    for cv_idx = 1:n_folds
        test_set = find (dt_inchi_mf_fp.cv == cv_idx);
        train_set = find (dt_inchi_mf_fp.cv ~= cv_idx);
        
        KX_list_train_test = cellfun(@(x) x(train_set, test_set), KX_list, ...
            'UniformOutput', false);
        KX_list_test = cellfun(@(x) x(test_set, test_set), KX_list, ...
            'UniformOutput', false);
        Y_train = Y(: ,train_set);
        Y_C_test = Y_C.getSubset (test_set);
        
        load ([outputDir, 'model_', MP_IOKR_Defaults.param2str(param), ...
            '_cv=', num2str(cv_idx), '.mat'], 'train_model');
        scores(test_set) = Test_IOKR (KX_list_train_test, KX_list_test, train_model, ...
            Y_train, Y_C_test, train_model.center);    
    end % for
    
    %--------------------------------------------------------------
    % Computation of the ranks
    %--------------------------------------------------------------
    inchi = dt_inchi_mf_fp.inchi_key_1;
    ranks = getRanksBasedOnScores (Y_C, inchi, scores);
    rank_perc = getRankPerc (ranks, arrayfun (@(id) Y_C.getCandidateSet (id, false, 'num'), ...
        1:Y_C.getNumberOfExamples()));
    rank_perc_100 = rank_perc(1:100);
    
    disp (round (rank_perc_100([1, 5, 10, 20]), 3));

    %--------------------------------------------------------------
    % Write out performance for comparison (top-k accuracy, ranks)
    %--------------------------------------------------------------
    cmp = dt_inchi_mf_fp.mol_id;
    ranks_t = table (cmp, ranks);
    filename = [outputDir 'ranks_cv_', MP_IOKR_Defaults.param2str(param), '.csv'];
    writetable (ranks_t, filename);
    
    filename = [outputDir 'rank_perc_100_cv', MP_IOKR_Defaults.param2str(param), '.txt'];
    save (filename,'rank_perc_100','-ascii');    
end % function 
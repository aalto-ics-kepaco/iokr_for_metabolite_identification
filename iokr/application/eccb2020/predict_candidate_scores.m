function [ ] = predict_candidate_scores(inputDir, independentSet, fps_def, param, outputDir)

    % ---------------------------------------------------
    % Handle default parameters
    % ---------------------------------------------------
    if nargin < 4
        param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
            {'debug_param', 'opt_param', 'iokr_param', 'data_param', 'ky_param'});
    
        % Number of cross-validation folds for the hyper-parameter
        % optimization, e.g. finding lambda (regularization parameter).

        % Output kernel: kappa_y
        param.ky_param.representation  = 'kernel';
        param.ky_param.type            = 'gaussian';
        param.ky_param.base_kernel     = 'tanimoto';
        param.ky_param.param_selection = 'entropy';

        param.iokr_param.model_representation = 'Chol_decomp_of_C';
    else
        param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, ...
            {'debug_param', 'opt_param', 'iokr_param', 'data_param', 'ky_param'});
    end % if 
    
    if nargin < 5
        outputDir = inputDir;
    end % if
    
    % ---------------------------------------------------
    % Load information IOKR model
    % ---------------------------------------------------
    model_fn = fullfile(inputDir, 'iokr_models', strcat(MP_IOKR_Defaults.param2str(param), '.mat')); 
    if ~ exist(model_fn, 'file')
        error('predict_candidate_scores:NoSuchFileOrDirectory', ...
            'Cannot find pre-trainined IOKR model: "%s"', model_fn);
    end % if
    tmp = load(model_fn);
    iokr_model = tmp.model;
    n_fps = tmp.n_fps;
    Y_train = tmp.Y_train;
    param = tmp.param;
    
    % ---------------------------------------------------
    % Load information of the independent dataset
    % ---------------------------------------------------
    indp_dir = fullfile(inputDir, 'independent', independentSet);
    cmps_fn = fullfile(indp_dir, 'compounds.csv'); 
    if ~ exist(cmps_fn, 'file')
        error('predict_candidate_scores:NoSuchFileOrDirectory', ...
            'Cannot find the compound list for the independent dataset: "%s". Use "fingerID list-independent-compounds"', ...
            cmps_fn);
    end % if
    cmps = readtable(cmps_fn, 'Delimiter', '\t', 'ReadVariableNames', false);
        
    % ---------------------------------------------------
    % Check availability of the candidate DB
    % ---------------------------------------------------
    % TODO: That should be solved more flexible.
    db_fn = fullfile(inputDir, 'independent', independentSet, 'candidates.db');
    if ~ exist(db_fn, 'file')
        error('predict_candidate_scores:NoSuchFileOrDirectory', ...
            'Cannot find the candidate DB: "%s"', db_fn);
    end % if
    db = sqlite(db_fn, 'readonly');
    if ~ db.IsOpen
        error('predict_candidate_scores:DBError', 'Cannot open connection to DB: "%s"', db_fn);
    end % if
    
    % ---------------------------------------------------
    % Prepare output directory for the scores
    % ---------------------------------------------------
    outputDir = fullfile(outputDir, 'independent', independentSet, 'scores', ...
        MP_IOKR_Defaults.param2str(param));
    if ~ exist(outputDir)
        mkdir(outputDir);
    end % if
    
    % ---------------------------------------------------
    % Predict the candidate scores
    % ---------------------------------------------------
    for idx = 1:5  %size(cmps)
        spec_id = cmps.Var1{idx};
        fprintf('Process "%s" (%d/%d)\n', spec_id, idx, size(cmps, 1))
        
        % Load the train-test and test-test kernels
        % ---------------------------------------------------
        
        % NOTE: We are only dealing with the most highest scores
        %       fragmentation tree here.
        train_test_pattern = fullfile(indp_dir, 'kernels', sprintf('%s__01__*.train_test_kernels', spec_id));
        train_test_fn = dir(train_test_pattern);
        if isempty(train_test_fn)
            error('predict_candidate_scores:NoSuchFileOrDirectory', ...
                'Cannot find train-test-kernel file for spectra: "%s"', spec_id);
        end % if
        [KX_trts, ~] = read_challenge_kernel_( ...
            fullfile(indp_dir, 'kernels', train_test_fn.name), ...
            param.data_param.availInputKernels);
        
        test_test_pattern = fullfile(indp_dir, 'kernels', sprintf('%s__01__*.test_test_kernels', spec_id));
        test_test_fn = dir(test_test_pattern);
        if isempty(test_test_fn)
            error('predict_candidate_scores:NoSuchFileOrDirectory', ...
                'Cannot find test-test-kernel file for spectra: "%s"', spec_id);
        end % if
        [KX_tsts, ~] = read_challenge_kernel_( ...
            fullfile(indp_dir, 'kernels', test_test_fn.name), ...
            param.data_param.availInputKernels);
        
        % Load the candidate information and fingerprints from the database
        % -----------------------------------------------------------------
        sqlstmt = strcat("SELECT molecule, %s FROM fingerprints_data", ...
                         "   INNER JOIN candidates_spectra cs ", ...
                         "      ON fingerprints_data.molecule = cs.candidate", ...
                         "   WHERE spectrum IS '%s'");
        rows = fetch(db, sprintf(sqlstmt, fps_def, spec_id));
        n_cand = size(rows, 1);
        fprintf('Number of candidates: %d\n', n_cand);
        
        cnd_ids = cell(n_cand, 1);
        Y_Ci = zeros(n_fps, n_cand);
        for j = 1:n_cand
            cnd_ids{j} = rows{j, 1};
            Y_Ci(fps_str_to_vec(rows{j, 2}) + 1, j) = 1;
        end % for
        Y_Ci = CandidateSets(DataHandle(struct('data', Y_Ci, 'id', {cnd_ids'}, 'num', n_cand)), 1);
        
        % Predict the scores using the pre-trained IOKR model
        % ---------------------------------------------------
        [scores_test, ~, hh_test] = Test_IOKR (KX_trts, KX_tsts, iokr_model, ...
                Y_train, Y_Ci, param.iokr_param.center);
        
        % Write out the scores    
        % --------------------
        res = table(cnd_ids, scores_test{1}', scores_test{1}' / hh_test, repmat (hh_test, [numel(scores_test{1}), 1]), ...
            'VariableNames', {'inchi', 'score', 'normalized_score', 'norm'});
        writetable(res, fullfile(outputDir, strcat(spec_id, '.csv')), ...
            'Delimiter', '\t', 'QuoteStrings', false);
        
        % Print out some statistics about the correct structures
        % ------------------------------------------------------
        rows = fetch(db, sprintf("SELECT molecule FROM spectra WHERE name is '%s'", spec_id));
        ranks = tiedrank(- scores_test{1}');
        disp(ranks);
        correct_idx = Y_Ci.findExampleInCandidateSet(1, rows{1, 1});
        fprintf('Rank of the correct structure (average): %d\n', ranks(correct_idx));
    end % for
        
    close(db);
end % function

function fps_vec = fps_str_to_vec(fps_str)
    fps_str_splitted = strsplit(fps_str, ',');
    fps_vec = cellfun (@(x) str2double(x), fps_str_splitted);
end % function

function [KX, KX_names] = read_challenge_kernel_(filename, input_kernel_names, split_char, verbose)
    %% READ_CHALLENGE_KERNEL reads in kernels line by line
    % Function to read training-vs-test or test vs test kernels calculated
    % with the CSI:FingerID framework.
    %
    % Structure of the training-vs-test kernel files:
    %   KERNEL_1[split_char]k_1(x_1, x')[split_char]k_1(x_2, x')...[split_char]k_1(x_n, x')\n
    %   KERNEL_2[split_char]k_2(x_1, x')[split_char]k_2(x_2, x')...[split_char]k_2(x_n, x')\n
    %   ...
    %
    %   with n being the number of training examples and x' being the
    %   specific test example.
    %
    %   E.g.
    %   PPKr 0.12 0.23 ... 0.31
    %   MLIP 0.52 0.4233 ... 0.769
    %   ...
    if nargin < 3
        split_char = ' ';
    end % if
    
    if nargin < 4
        verbose = false;
    end % if

    if verbose
        tic;
    end % if

    fid = fopen(filename, 'r');
    if (fid < 0)
        error('read_challenge_kernel_:NoSuchFileOrDirectory', 'Cannot open file: %s', filename);
    end % if
    
    input_kernel_names = upper(input_kernel_names);
    KX_names = cell(0);
    KX = cell(0);
    k = 1;
    
    line = fgetl(fid);
    while (ischar(line))
        line_splitted = strsplit(line, split_char);
        
        if ismember(upper(line_splitted{1}), input_kernel_names)
            KX_names{k} = upper(line_splitted{1});
            KX{k} = cellfun (@(x) str2double(x), line_splitted(2:end))';
            k = k + 1;
        end % if
        
        line = fgetl(fid);
        
    end % if
    
    fclose(fid);
    
    % Need to ensure, that the order of the kernels matches the order used during training of the IOKR model     
    [~, locb] = ismember(input_kernel_names, KX_names);
    assert (all(locb > 0), 'Could not find all kernel names!')
    KX_names = KX_names(locb);
    assert (all(cellfun(@strcmp, KX_names, input_kernel_names')))
    KX = KX(locb);
    
    if verbose
        toc;
    end % if
end % function 

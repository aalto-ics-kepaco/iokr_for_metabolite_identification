function run_IOKR_independent_set(base_dir_training, base_dir_indep)
%% RUN_IOKR_INDEPENDENT_SET Scoring script for independent set of spectra
% Requires preprocessed data: path/to/independent_set
%   - kernels/      : kernel files between independent MSMS and training
%                     data
%   - candidates/   : csv-files containing candidates, fingerprints for all
%                     candidates

    %% Set up parameters
    %% -----------------
    param = MP_IOKR_Defaults.setDefaultsIfNeeded(struct(), {'data_param','ky_param','iokr_param', 'opt_param'});
    param.iokr_param.model_representation = 'only_C'; % 'Chol_decomp_of_C'
    param.indep_param = struct('remove_test_inchikey2D_from_training', true);
    param.data_param.inputKernel = '__ALL__';
    param.data_param.use_fps_mask = true;
    
    %% Set up directory paths
    %% ----------------------
    kernel_dir = strcat(base_dir_training, '/kernels/');
    raw_fps_dir = strcat(base_dir_training, '/fingerprints/');
    
    mat_fps_dir = strcat(base_dir_indep, '/fingerprints/');
    if ~ (mkdir(mat_fps_dir) == 1)
        error('run_IOKR_independent_set:CanNotCreateDirectory', 'Could not create the independent set fps directory "%s".', ...
            mat_fps_dir);
    end % if
    
    %% Define output directory
    %% -----------------------
    % find available kernels:Check .mat-files first than .txt-files
    kernel_files_ext = '.mat';
    kernel_files = dir(strcat(kernel_dir, '/*', kernel_files_ext));
    if isempty(kernel_files)
        kernel_files_ext = '.txt';
        kernel_files = dir(strcat(kernel_dir, '/*', kernel_files_ext));
        
        if isempty(kernel_files)
            error('run_IOKR_independent_set:NoKernelsAvailable', ...
                'Directory "%s" does not contain kernel files.', kernel_dir);
        end % if
    end % if
    
    % get names of available input kernels
    param.data_param.availInputKernels = ... 
        arrayfun(@(file) basename(file.name), kernel_files, 'UniformOutput', false);
    
    % hash model parameters
    model_param_hash = string2hash( ...
        strcat(MP_IOKR_Defaults.param2str(param, true, {'ky_param', 'iokr_param','data_param'}), ...
            '_rm-test-inchikeys=', num2str(param.indep_param.remove_test_inchikey2D_from_training), ...
            '_use-fps-mask', num2str(param.data_param.use_fps_mask)));
    output_dir = sprintf('%s/msms_scores_%.0f/', base_dir_indep, model_param_hash);
    if ~ (mkdir(output_dir) == 1)
        error('run_IOKR_independent_set:CannotCreateDirectory', 'Could not create the output directory "%s".', ...
            output_dir);
    end % if
    
    % write out model parameter as json-file
    model_param_json_fn = sprintf('%s/model_parameters.json', output_dir);
    fid = fopen(model_param_json_fn, 'w');
    if fid == -1
        error('run_IOKR_independent_set:CannotOpenFile', 'Cannot open model parameter json-file "%s" for writing.', ...
            model_param_json_fn);
    end % if
    fprintf(fid, '%s', jsonencode(param));
    fclose(fid);
    
    %% Check whether trained model is available
    %% ----------------------------------------
    % Load the list of independent compounds. This list can be created
    % using the csi-fingerid CLI: 'fingerID list-independent-compounds --indep-set=INDEPSET > independent/INDEPSET/compounds'
    % NOTE: If the compounds are not known, i.e. no InChI or SMILES is
    %       available in the ms-file. Than the table contains only the
    %       'spec_id', e.g. the 'inchikey2D' column would be filled with
    %       '__NOT__KNOWN__'.
    cmps_indep = loadCompoundList(base_dir_indep, 'csi_fingerid');
    
    % Load the list of training compounds. This list can be created using
    % the csi-fingerid CLI: 'fingerID list-compounds > compounds'
    cmps_train = loadCompoundList(base_dir_training, 'csi_fingerid');

    model_fn = strcat(output_dir, '/iokr_model.mat');
    if ~ exist(model_fn, 'file')
        fprintf('Train model ...\n');
        
        %% Remove compounds from training that are in the independent set
        if param.indep_param.remove_test_inchikey2D_from_training        
            % Exclude spectra from training with molecular structures, based on
            % the 2D inchikey, present in the training set.
            cmp_is_used_for_training = ~ ismember(cmps_train.inchikey2D, cmps_indep.inchikey2D);
        else
            % Use all compounds
            cmp_is_used_for_training = true(size(cmps_train, 1), 1);
        end % if
        
        %% Load training data
        % load kernels 
        KX_list_train = loadInputKernelsIntoList(kernel_dir, param, kernel_files_ext);
        if strcmp(kernel_files_ext, '.txt')
            % store kernel files in .mat-format
            for idx = 1:length(KX_list_train)
                saveKernel(strcat(kernel_dir, '/', param.data_param.availInputKernels{idx}, '.mat'), ...
                    KX_list_train{idx}, cmps_train.spec_id);
            end % if
        end % if    
        KX_list_train = cellfun(@(KX) KX(cmp_is_used_for_training, cmp_is_used_for_training), KX_list_train, ...
            'UniformOutput', false);

        % Load training fingerprints
        % NOTE: Y is a (d x n)-matrix, with d being the dimension of the fps,
        %       and n being the number of examples.
        % NOTE: Only those fingerprints are loaded that are needed for the
        %       model.
        [Y_train, fps_mask] = loadFingerprints(raw_fps_dir, mat_fps_dir, cmps_train.spec_id(cmp_is_used_for_training));
        if param.data_param.use_fps_mask
            Y_train = Y_train(fps_mask, :); 
        end % if
   
        %% Train IOKR
        iokr_model = Train_IOKR(KX_list_train, Y_train, param.ky_param, param.opt_param, param.iokr_param, ...
            true);
        
        iokr_model.cmp_is_used_for_training = cmp_is_used_for_training;
        
        save(model_fn, 'iokr_model', '-v7.3');
    else
        fprintf('Load model ...\n');
        
        load(model_fn, 'iokr_model');
        
        cmp_is_used_for_training = iokr_model.cmp_is_used_for_training;
        
        % Load training fingerprints
        % NOTE: Y is a (d x n)-matrix, with d being the dimension of the fps,
        %       and n being the number of examples.
        % NOTE: Only those fingerprints are loaded that are needed for the
        %       model.
        [Y_train, fps_mask] = loadFingerprints(raw_fps_dir, mat_fps_dir, cmps_train.spec_id(cmp_is_used_for_training));
        if param.data_param.use_fps_mask
            Y_train = Y_train(fps_mask, :); 
        end % if
    end % if


    %% Apply IOKR model: Predict candidate scores for the different candidate sets
    %% ---------------------------------------------------------------------------
    
    %% 1) Create candidate set object
    indep_cand_lists = cell(size(cmps_indep, 1), 1);
    for idx = 1:length(indep_cand_lists)
        tmp = dir(strcat(base_dir_indep, '/candidates/', cmps_indep.spec_id{idx}, '__*__*.candidates'));
        if isempty(tmp)
            warning('Cannot find candidate set for "%s".', cmps_indep.spec_id{idx});
        elseif length(tmp) > 1
            error('run_IOKR_independent_set:NotImplementedYet', ...
                'Currently only one candidate set per spectra can be processed.');
        end % if
        
        indep_cand_lists{idx} = basename(tmp(1).name);
    end % for
    clear tmp;
    
    if param.data_param.use_fps_mask
        Y_C = CandidateSetsFile(strcat(base_dir_indep, '/candidates/'), indep_cand_lists, fps_mask);
    else
        Y_C = CandidateSetsFile(strcat(base_dir_indep, '/candidates/'), indep_cand_lists, '__ALL__');
    end % if
    
    %% 2) Load the test kernels
    KX_list_train_test = arrayfun(@(c) nan(sum(cmp_is_used_for_training), size(cmps_indep, 1)), ... 
        1:length(param.data_param.availInputKernels), 'UniformOutput', false);
    KX_list_test_test = arrayfun(@(c) zeros(size(cmps_indep, 1)), ... 
        1:length(param.data_param.availInputKernels), 'UniformOutput', false);
    
    kernel_dir_indep = strcat(base_dir_indep, '/kernels/');
    for idx = 1:length(indep_cand_lists)        
        % Load kernels
        [tmp_KX_train_test, KX_names_train_test] = read_challenge_kernel( ...
            strcat(kernel_dir_indep, basename(indep_cand_lists{idx}), '.train_test_kernels'), ' ', true);
        [tmp_KX_test_test, KX_names_test_test] = read_challenge_kernel( ...
            strcat(kernel_dir_indep, basename(indep_cand_lists{idx}), '.test_test_kernels'), ' ', true);

        % Check order of the test kernels against the ones used for
        % the model training.
        [~, locb] = ismember(upper(param.data_param.availInputKernels), upper(KX_names_train_test));
        if ~issorted(locb)
            warning('Train-vs-test: "locb" is not sorted.')
        end % if
        tmp_KX_train_test = tmp_KX_train_test(locb);
        tmp_KX_train_test = cellfun(@(KX) KX(cmp_is_used_for_training, :), tmp_KX_train_test, ...
            'UniformOutput', false);

        [~, locb] = ismember(upper(param.data_param.availInputKernels), upper(KX_names_test_test));
        if ~issorted(locb)
            warning('Test-vs-test: "locb" is not sorted.')
        end % if
        tmp_KX_test_test = tmp_KX_test_test(locb);
        
        for k_idx = 1:length(param.data_param.availInputKernels)
            KX_list_train_test{k_idx}(:, idx) = tmp_KX_train_test{k_idx};
            KX_list_test_test{k_idx}(idx, idx) = tmp_KX_test_test{k_idx};
        end % for
    end % for
    
    %% 3) Predict candidate scores 
    [scores, ~, squared_norm_of_h] = Test_IOKR(KX_list_train_test, KX_list_test_test, iokr_model, Y_train, ...
        Y_C, iokr_model.ker_center);
    
    %% 4) Store results
    ordinal_ranks_true_candidates = zeros(length(indep_cand_lists), 1);
    average_ranks_true_candidates = zeros(length(indep_cand_lists), 1);
    dense_ranks_true_candidates = zeros(length(indep_cand_lists), 1);
    min_ranks_true_candidates = zeros(length(indep_cand_lists), 1);
    max_ranks_true_candidates = zeros(length(indep_cand_lists), 1);
    
    for idx = 1:length(indep_cand_lists)
        inchi2D = Y_C.getCandidateSet(idx, false, 'inchi2D');
        inchikey2D = Y_C.getCandidateSet(idx, false, 'id');
        is_true_candidate = Y_C.findExampleInCandidateSet(idx, cmps_indep.inchikey2D{idx});
        
        % "Random rank": two molecules with the same score, more or less
        % randomly ordered. This is because the order of the candidate in
        % the set is not sorted.
        %
        % PYTHON: scipy.stats.rankdata(score, 'ordinal')
        [~, tmp] = sort(scores{idx}, 'descend');
        ordinal_ranks_based_on_msms_score = zeros(length(inchi2D), 1);
        ordinal_ranks_based_on_msms_score(tmp) = 1:length(inchi2D);
        
        ordinal_ranks_true_candidates(idx) = ordinal_ranks_based_on_msms_score(is_true_candidate);
        
        % "Mean rank": two molecules with the same score will have
        % averaged rank. 
        %
        % PYTHON: scipy.stats.rankdata(score, 'average')
        average_ranks_based_on_msms_score = tiedrank(-scores{idx})';
        
        average_ranks_true_candidates(idx) = average_ranks_based_on_msms_score(is_true_candidate);
        
        % "Shared rank (dense)": two molecules with the same score will have the
        % same rank.  
        %
        % PYTHON: scipy.stats.rankdata(score, 'dense')
        [~, ~, dense_shared_ranks_based_on_msms_score] = unique(-scores{idx});
        
        dense_ranks_true_candidates(idx) = dense_shared_ranks_based_on_msms_score(is_true_candidate);
        
        % "Shared rank (gaps, min)": two molecules with the same score will have the
        % same rank.  
        %
        % PYTHON: scipy.stats.rankdata(score, 'min')
        min_shared_ranks_based_on_msms_score = floor(tiedrank(-scores{idx}))';
        
        min_ranks_true_candidates(idx) = min_shared_ranks_based_on_msms_score(is_true_candidate);
        
        % "Shared rank (gaps, max)": two molecules with the same score will have the
        % same rank.  
        %
        % PYTHON: scipy.stats.rankdata(score, 'max')
        max_shared_ranks_based_on_msms_score = ceil(tiedrank(-scores{idx}))';
        
        max_ranks_true_candidates(idx) = max_shared_ranks_based_on_msms_score(is_true_candidate);
        
        % Set up score talbe for the current candidate set
        score_table = table( ...
            inchi2D, ...
            inchikey2D, ...
            scores{idx}', ...
            repmat(squared_norm_of_h(idx), [length(scores{idx}), 1]), ...
            is_true_candidate, ...
            ordinal_ranks_based_on_msms_score, ...
            average_ranks_based_on_msms_score, ...
            dense_shared_ranks_based_on_msms_score, ...
            min_shared_ranks_based_on_msms_score, ...
            max_shared_ranks_based_on_msms_score, ...
            'VariableNames', ...
            {'inchi2D', ...
            'inchikey2D', ...
            'msms_score', ...
            'squared_norm_of_hx', ...
            'is_true_candidate', ...
            'ordinal_rank_based_on_msms_score', ...
            'average_rank_based_on_msms_score', ...
            'dense_rank_based_on_msms_score', ...
            'min_rank_based_on_msms_score', ...
            'max_rank_based_on_msms_score'});
        
        score_table_fn = sprintf('%s/%s.csv', output_dir, basename(indep_cand_lists{idx}));
        writetable(score_table, score_table_fn, 'Delimiter', '\t')
    end % for
    
    %% Determine rank percentages
    max_cand_num = max(arrayfun(@(idx) Y_C.getCandidateSet(idx, false, 'num'), 1:Y_C.getNumberOfExamples()));
    ordinal_rank_perc = getRankPerc(ordinal_ranks_true_candidates, max_cand_num);
    average_rank_perc = getRankPerc(average_ranks_true_candidates, max_cand_num);
    dense_rank_perc = getRankPerc(dense_ranks_true_candidates, max_cand_num);
    min_rank_perc = getRankPerc(min_ranks_true_candidates, max_cand_num);
    max_rank_perc = getRankPerc(max_ranks_true_candidates, max_cand_num);
    
    topk = 1:100;
    rank_perc_table = table(topk', ordinal_rank_perc(topk), average_rank_perc(topk), dense_rank_perc(topk), ...
        min_rank_perc(topk), max_rank_perc(topk), 'VariableNames', {'topk', 'ordinal', 'average', 'dense', ...
        'min', 'max'});
    
    rank_perc_table_fn = sprintf('%s/rank_perc.csv', output_dir);
    writetable(rank_perc_table, rank_perc_table_fn, 'Delimiter', '\t')
end % function


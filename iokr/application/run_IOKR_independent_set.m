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
    param.iokr_param.model_representation = 'Chol_decomp_of_C';
    
    param.indep_param = struct('remove_test_inchikey2D_from_training', true);
    
    %% Load training data
    %% ------------------

    % Load the list of training compounds. This list can be created using
    % the csi-fingerid CLI: 'fingerID list-compounds > compounds'
    cmps_train = loadCompoundList(base_dir_training, 'csi_fingerid');
    
    % Load training kernels
    kernel_dir = strcat(base_dir_training, '/kernels/');
    
    % try to load .mat-file kernels
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
    param.data_param.inputKernel = '__ALL__';
    
    % load kernels 
    KX_list_train = loadInputKernelsIntoList(kernel_dir, param, kernel_files_ext);
    if strcmp(kernel_files_ext, '.txt')
        % store kernel files in .mat-format
        for idx = 1:length(KX_list_train)
            saveKernel(strcat(kernel_dir, '/', param.data_param.availInputKernels{idx}, '.mat'), ...
                KX_list_train{idx}, cmps_train.spec_id);
        end % if
    end % if
     
    % Load training fingerprints
    % NOTE: Y is a (d x n)-matrix, with d being the dimension of the fps,
    %       and n being the number of examples.
    tmp = loadFingerprints(strcat(base_dir_training, '/fingerprints/'), cmps_train.spec_id);
    Y_train = tmp.fps;
    fps_mask = tmp.mask; % TODO: What should we do with the fingerprint mask?
    clear tmp;
    
    %% Define output directory
    %% -----------------------
    model_param_hash = string2hash( ...
        strcat(MP_IOKR_Defaults.param2str(param, true, {'ky_param', 'iokr_param','data_param'}), ...
            '_rem-test-inchikeys=', num2str(param.indep_param.remove_test_inchikey2D_from_training)), ...
            'sdbm');
    output_dir = sprintf('%s/msms_scores_%.0f/', base_dir_indep, model_param_hash);
    if ~ (mkdir(output_dir) == 1)
        error('run_IOKR_independent_set:CanNotCreateDirectory', 'Could not create the output directory "%s".', ...
            output_dir);
    end % if
    
    %% Remove compounds from training that are in the independent set
    %% --------------------------------------------------------------
    if param.indep_param.remove_test_inchikey2D_from_training        
        % Load the list of independent compounds. This list can be created
        % using the csi-fingerid CLI: 'fingerID list-independent-compounds --indep-set=INDEPSET > independent/INDEPSET/compounds'
        % NOTE: If the compounds are not known, i.e. no InChI or SMILES is
        %       available in the ms-file. Than the table contains only the
        %       'spec_id', e.g. the 'inchikey2D' column would be filled with
        %       '__NOT__KNOWN__'.
        cmps_indep = loadCompoundList(base_dir_indep, 'csi_fingerid');
        
        % Exclude spectra from training with molecular structures, based on
        % the 2D inchikey, present in the training set.
        cmp_is_used_for_training = ~ ismember(cmps_train.inchikey2D, cmps_indep.inchikey2D);
    else
        % Use all compounds
        cmp_is_used_for_training = true(size(cmps_train, 1), 1);
    end % if
    
    KX_list_train = cellfun(@(KX) KX(cmp_is_used_for_training, cmp_is_used_for_training), KX_list_train, ...
        'UniformOutput', false);
    Y_train = Y_train(:, cmp_is_used_for_training);
    
    %% Train IOKR
    %% ----------
    model_fn = strcat(output_dir, '/iokr_model.mat');
    if ~ exist(model_fn, 'file')
        iokr_model = Train_IOKR(KX_list_train, Y_train, param.ky_param, param.opt_param, param.iokr_param);
        save(model_fn, 'iokr_model', '-v7.3');
    else
        load(model_fn, 'iokr_model');
        iokr_model = iokr_model.iokr_model;
    end % if

    %% Apply IOKR model: Predict candidate scores for the different candidate sets
    %% ---------------------------------------------------------------------------
    kernel_dir_indep = strcat(base_dir_indep, '/kernels/');
        
    % Create candidate set object
%     lut = arrayfun(@(x) basename(x.name), dir(strcat(base_dir_indep, '/candidates/*.candidates')), 'UniformOutput', false);
%     Y_C = CandidateSetsFile(strcat(base_dir_indep, '/candidates/'), lut, fps_mask);
    
    for indep_spec_id = cmps_indep.spec_id'
        indep_cand_lists = dir(strcat(base_dir_indep, '/candidates/', indep_spec_id{1}, '__*__*.candidates'));
        
        fprintf('Spectra ID: "%s", with %d candidate list(s).\n', indep_spec_id{1}, length(indep_cand_lists))
        
        for indep_cand_list = indep_cand_lists'
            % Load kernels
            [KX_list_train_test, KX_names_train_test] = read_challenge_kernel( ...
                strcat(kernel_dir_indep, basename(indep_cand_list.name), '.train_test_kernels'), ' ', true);
            [KX_list_test_test, KX_names_test_test] = read_challenge_kernel( ...
                strcat(kernel_dir_indep, basename(indep_cand_list.name), '.test_test_kernels'), ' ', true);
            
            % Check order of the test kernels against the ones used for
            % the model training.
            [~, locb] = ismember(upper(param.data_param.availInputKernels), upper(KX_names_train_test));
            if ~issorted(locb)
                warning('Train-vs-test: "locb" is not sorted.')
            end % if
            KX_list_train_test = KX_list_train_test(locb);
            KX_list_train_test = cellfun(@(KX) KX(cmp_is_used_for_training, :), KX_list_train_test, ...
                'UniformOutput', false);
            
            [~, locb] = ismember(upper(param.data_param.availInputKernels), upper(KX_names_test_test));
            if ~issorted(locb)
                warning('Test-vs-test: "locb" is not sorted.')
            end % if
            KX_list_test_test = KX_list_test_test(locb);
            
            % Predict candidate scores 
            Y_C = CandidateSetsFile(strcat(base_dir_indep, '/candidates/'), {basename(indep_cand_list.name)}, '__ALL__');
            [scores, ~, squared_norm_of_h] = Test_IOKR(KX_list_train_test, KX_list_test_test, iokr_model, Y_train, ...
                Y_C, iokr_model.ker_center);
            
            assert(~ any(isnan(scores{1})));
            assert(length(scores) == 1);
            
            
        end % for
    end % for
end % function


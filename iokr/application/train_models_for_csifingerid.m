function train_models_for_csifingerid (inputDir, outputDir)
%======================================================
% DESCRIPTION:
% Script for running IOKR
%
% INPUTS:
% inputDir:     directory in which the data are contained
% result_dir:   directory in which the results will be saved
%
%======================================================
    if (~ exist (outputDir, 'dir'))
        mkdir (outputDir);
    end % if

    %--------------------------------------------------------------
    % Load Data
    %--------------------------------------------------------------
    % inchi keys, molecular formulas, fingerprints
    load ([inputDir '/compound_info.mat'], 'dt_inchi_mf_fp'); 
    Y     = full (dt_inchi_mf_fp.fp_full)';
    
    % Fingerprint mask 
    fp_mask = load_fingerprint_mask ([inputDir, '/fingerprints.mask']);
    Y = Y(fp_mask, :);
    
    % Input kernels
    kernel_files = dir ([inputDir '/kernels/*.txt']);

    % Parameters
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
        {'opt_param', 'iokr_param', 'data_param', 'ky_param'});       
    param.data_param.availInputKernels = arrayfun (@(file) basename (file.name), ...
        kernel_files, 'UniformOutput', false);
    
    % Convert kernel matrices to matfiles, if those do not exist yet.
    if isempty (dir ([inputDir, '/kernels/*.mat']))
        convert_kernel_files_to_mat ([inputDir, '/kernels/']);
    end % if 
    
    KX_list = loadInputKernelsIntoList ([inputDir '/kernels/'], param, '.mat');
    
    %--------------------------------------------------------------
    % Model training using cross-validation
    %--------------------------------------------------------------
    disp ('Model training')
    
    filename = [outputDir, 'model_', MP_IOKR_Defaults.param2str(param), '.mat'];
    if (exist (filename, 'file'))
        load (filename, 'train_model');
    else
        train_model = Train_IOKR (KX_list, Y, ...
        param.ky_param, param.opt_param, param.iokr_param);
    
        switch train_model.KY_par.representation
            case 'feature'
                [~, process_output] = output_feature_preprocessing_train (Y, ...
                    train_model.ker_center);
            case 'kernel'
                [~, process_output] = output_kernel_preprocessing_train (Y, ...
                    train_model.KY_par, train_model.ker_center);
        end % switch 
        train_model.process_output = process_output;
        train_model.kernel_names = param.data_param.availInputKernels;
        
        save (filename, 'train_model', '-v7.3');
    end % if
    
    %--------------------------------------------------------------
    % Save the model: As binary for csi-fingerid 
    %--------------------------------------------------------------
    if (~ exist ([outputDir, '/model_as_binary/'], 'dir'))
        mkdir ([outputDir, '/model_as_binary/']);
    end % if
    
    write_out_iokr_model_as_binary_files (train_model, [outputDir, '/model_as_binary/']);
    
    %--------------------------------------------------------------
    % Prediction & scoring 
    %--------------------------------------------------------------
    error ('Not implemented yet.')
    
    disp ('Prediction & Scoring')
    
    mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, cand);
    load ([inputDir, '/cand.mat'], 'cand');
    
    Y_C = CandidateSets (DataHandle (cand), mf_corres, 'ALL', fp_mask);
    scores = Test_IOKR (KX_list, KX_list, train_model, Y, Y_C, train_model.center);
    
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
    filename = [outputDir 'ranks_', MP_IOKR_Defaults.param2str(param), '.csv'];
    writetable (ranks_t, filename);
    
    filename = [outputDir 'rank_perc_100', MP_IOKR_Defaults.param2str(param), '.txt'];
    save (filename,'rank_perc_100','-ascii');    
end


%     filename = [outputDir 'scores_reclass_mkl=' iokr_param.mkl ...
%         '_kernel=' ky_param.type ...
%         '_base=' ky_param.base_kernel '_' ky_param.param_selection ...
%         '_model_representation=', iokr_param.model_representation, '.mat'];
%     save (filename,'scores','-v7.3');
% 
%     
%     rank = NaN(n,1);
%     cand_num = zeros(n,1); % vector containing the number of candidates for each test example
%     for j = 1:n
%         inchi_c = cand(mf_corres(j)).id;
%         [~,IX] = sort(score{j},'descend');
% 
%         rank(j) = find(strcmp(inchi_c(IX), inchi{j}));
%         cand_num(j) = length(score{j});
%     end
%         
%     % Computation of the percentage of identified metabolites in the top k
%     nel = hist(rank, 1:max(cand_num));
%     rank_perc = cumsum(nel)';
%     rank_perc = rank_perc/n*100;
%     rank_perc_100 = rank_perc(1:100);

% filename = [outputDir '/model_reclass_mkl=' iokr_param.mkl ...
%         '_kernel=' ky_param.type ...
%         '_base=' ky_param.base_kernel '_' ky_param.param_selection ...
%         '_model_representation=', iokr_param.model_representation, '.mat'];
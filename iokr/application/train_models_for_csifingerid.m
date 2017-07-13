function [rank, score, train_model] = train_models_for_csifingerid( ...
    inputDir, outputDir, cand )
%======================================================
% DESCRIPTION:
% Script for running IOKR
%
% INPUTS:
% inputDir:     directory in which the data are contained
% result_dir:   directory in which the results will be saved
%
%======================================================

    rng (10);
%     addpath(genpath('..'));

    %--------------------------------------------------------------
    % Load Data
    %--------------------------------------------------------------

    % inchi keys, molecular formulas, fingerprints
    load ([inputDir '/compound_info.mat'], 'dt_inchi_mf_fp');
    
    inchi = dt_inchi_mf_fp.inchi_key_1; 
    
    % Extract fingerprints
    Y = full (dt_inchi_mf_fp.fp_masked)';
    [~,n] = size(Y);

    % Candidates description
    % load([inputDir '/cand.mat'],'cand');
    
    mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, cand);

    % Input kernels
    kernel_files = dir ([inputDir '/*.txt']);
    KX_list = arrayfun (@(file) loadKernel ([inputDir '/' file.name]), kernel_files, ...
        'UniformOutput', false);
    
    % Parameters
    iokr_param = struct('center',1,'mkl','alignf',...
        'model_representation','Chol_decomp_of_C');
    select_param = struct( ...
        'cv_type','loocv', ...
        'lambda', [1e-5 1e-4 1e-3 1e-2 5e-2 1e-1 5e-1 1 10 100]);
    ky_param = struct( ...
        'type','gaussian', ...
        'base_kernel','tanimoto',...
        'param_selection','entropy');
    output_param = struct('representation','kernel','kernel_param',ky_param);

    %--------------------------------------------------------------
    % Cross-validation
    %--------------------------------------------------------------

    rank = NaN(n,1);
    cand_num = zeros(n,1); % vector containing the number of candidates for each test example

    
    train_model = Train_IOKR(KX_list, Y, output_param, ...
            select_param, iokr_param);

    % Prediction 
    Y_C = arrayfun(@(x) x.data, cand(mf_corres(1:n)),'UniformOutput',false);

    disp ('TESTING')
    [score, process_output] = Test_IOKR (KX_list, KX_list, train_model, ...
        Y, Y_C, iokr_param.center);
    
    train_model.process_output = process_output;

    % Computation of the ranks
    for j = 1:n
        inchi_c = cand(mf_corres(j)).id;
        [~,IX] = sort(score{j},'descend');

        rank(j) = find(strcmp(inchi_c(IX), inchi{j}));
        cand_num(j) = length(score{j});
    end
        
    % Computation of the percentage of identified metabolites in the top k
    nel = hist(rank, 1:max(cand_num));
    rank_perc = cumsum(nel)';
    rank_perc = rank_perc/n*100;
    rank_perc_100 = rank_perc(1:100);
    
    disp (round (rank_perc_100([1, 5, 10, 20]), 3));

    filename = [outputDir 'score_reclass_mkl=' iokr_param.mkl ...
        '_kernel=' ky_param.type ...
        '_base=' ky_param.base_kernel '_' ky_param.param_selection ...
        '_model_representation=', iokr_param.model_representation, '.mat'];
    save(filename,'score','-v7.3');
    
    filename = [outputDir 'rank_reclass_mkl=' iokr_param.mkl ...
        '_kernel=' ky_param.type ...
        '_base=' ky_param.base_kernel '_' ky_param.param_selection ...
        '_model_representation=', iokr_param.model_representation];
    save(filename,'rank','-ascii');
    
    filename = [outputDir 'model_reclass_mkl=' iokr_param.mkl ...
        '_kernel=' ky_param.type ...
        '_base=' ky_param.base_kernel '_' ky_param.param_selection ...
        '_model_representation=', iokr_param.model_representation, '.mat'];
    save (filename, 'train_model', '-v7.3');
end
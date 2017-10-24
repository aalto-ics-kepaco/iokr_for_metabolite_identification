function run_IOKR_for_candidate_scoring (input_dir, output_dir)
%% RUN_IOKR_FOR_CANDIDATE_SCORING Return the scoring for a set of spectra.
%   This script is a preprocessing step for Zheyangs project. For a list
%   of compounds and their spectra we train an IOKR model and score the 
%   candidates. The compounds are part of an artificial metabolic pathway. 
%
%   INPUTS:
%       inputDir            Directory path containing the input files:
%           - kernels/*.mat: Input kernels as mat-files
%           - compound_info.mat:
%                            Fingerprints, InChIs etc. of the training compounds
%           - cand.mat:      Candidate sets for all training compounds.
%           - pathway_compounds.csv:
%                            CSV-file containing the information about the 
%                            compounds in the pathway. Those compounds, are
%                            excluded from the training set and a scoring
%                            is provided for all their candidates.
%       outputDir           Directory path to store the scorings for all
%                           compounds in the pathway.
%       param               Parameter structure for the IOKR algorithm

    %--------------------------------------------------------------
    % Check the input arguments and set defaults
    %--------------------------------------------------------------
    if (nargin < 2)
        error ('run_IOKR_for_candidate_scoring:InvalidInput', ...
            'Not enough input arguments.');
    end % if
    if (~ exist (input_dir, 'dir'))
        error ('run_IOKR_for_candidate_scoring:InvalidInput', '%s: No such directory.', ...
            input_dir);
    end % if
    if (~ exist (output_dir, 'dir'))
        error ('run_IOKR_for_candidate_scoring:InvalidInput', '%s: No such directory.', ...
            output_dir);
    end % if
    
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
        {'debug_param', 'opt_param', 'iokr_param', 'data_param', 'ky_param'});    
    
    %--------------------------------------------------------------
    % Load and prepare data
    %--------------------------------------------------------------
    % Fingerprints (use the masked fps)
    load ([input_dir, '/compound_info.mat'], 'dt_inchi_mf_fp');
    Y = full (dt_inchi_mf_fp.fp_masked)';
    
    % Load the fingerprint mask. This does not effect the performance, but
    % speeds up the computation.
    fp_mask = load_fingerprint_mask ([input_dir, '/fingerprints.mask']);
    
    % Load the pathway compounds:
    %   These compounds needs to be excluded from the training.
    pathway_cmps = readtable ([input_dir, '/pathway_compounds.csv']);
    pathway_cmps_spec_id_all = strsplit (strjoin(pathway_cmps.spec_id_neg_all, ';'), ';');
    pathway_cmps_spec_id_all = pathway_cmps_spec_id_all(cellfun(@(c) numel(c), pathway_cmps_spec_id_all) > 0);
    
    % All the spectra which are associated with the any (3D, 2D, unchg 2D)
    % structures of the pathway compounds are excluded from the training.
    train_set = find (~ ismember (dt_inchi_mf_fp.mol_id, pathway_cmps_spec_id_all));
   
    fprintf ('Number of training examples: %d\n', length (train_set));
    
    % Load the input kernels
    kernel_dir = [input_dir, '/kernels/'];
    kernel_files = dir ([kernel_dir, '/*.mat']);
    param.data_param.availInputKernels = arrayfun (@(file) basename (file.name), kernel_files, ...
        'UniformOutput', false);
    KX_list = loadInputKernelsIntoList (kernel_dir, param, '.mat');
    
    %--------------------------------------------------------------
    % Train model
    %--------------------------------------------------------------
    KX_list_train = cellfun(@(x) x(train_set, train_set), KX_list, 'UniformOutput', false);
    Y_train       = Y(:, train_set);
    model_fn      = [output_dir, '/', MP_IOKR_Defaults.param2str(param), '.mat'];
    
    fprintf ('Train model ... ');
    tic;
    if (exist (model_fn, 'file'))
        fprintf ('load model from file ... ');
        % Load the train model
        load (model_fn, 'train_model');
    else
        train_model = Train_IOKR (KX_list_train, Y_train, ...
            param.ky_param, param.opt_param, param.iokr_param);
        % Save the train model
        save (model_fn, 'train_model', 'param');
    end % if    
    fprintf ('%.3fs\n', toc);
    fprintf ('Optimal lambda: %e\n', train_model.lambda_opt);
    
    %--------------------------------------------------------------
    % Score candidates for each test example
    %--------------------------------------------------------------
    % Load the candidates
    fprintf ('Load candiates ... ');
    tic,
        cand = load ([input_dir, '/cand_prepared.mat']);
    fprintf ('%.3fs\n', toc);
    
    % Get the corresponding candidates set each spectra in the test set:
    %   For some pathway compounds we do have several spectra. We predict
    %   the scores for each candidate for each spectra, so that we can
    %   later for example average the scores belonging to different spectra
    %   measured from the same compound.
    n_mol_pathway = size (pathway_cmps, 1);
    
    ranks_table = table;
    ranks_table_no_spectra = table;
    
    for idx_mol = 1:n_mol_pathway
        fprintf ('Scoring: %d/%d\n', idx_mol, n_mol_pathway);
        
        pathway_cmps_spec_id_best = strsplit (pathway_cmps.spec_id_neg{idx_mol}, ';');
        n_spec_best = numel (pathway_cmps_spec_id_best);
        
        if strcmp (pathway_cmps_spec_id_best{1}, '')
            % No spectra available: Pseudo-scoring is created
            
            score_list_fn = sprintf ('%s/scoring_list=%s_spec=%s_id=%s.csv', ...
                output_dir, ...
                pathway_cmps.molecular_formula{idx_mol}, ...
                'NOSPECTRA', ...
                pathway_cmps.inchikey{idx_mol});
            fprintf ('\tScoring-file: %s\n', score_list_fn);
            
            id1 = pathway_cmps.inchikey_1(idx_mol);
            score = 1;
            score_list = table (id1, score);
            
            writetable (score_list, score_list_fn, 'Delimiter', ',', ...
                'WriteRowNames', false, 'WriteVariableNames', true')
            
            id       = pathway_cmps.inchikey(idx_mol);
            rank     = 1;
            spec_id  = {'NOSPECTRA'};
            cand_num = -1;
            ranks_table_no_spectra = [ranks_table_no_spectra; table(id, rank, spec_id, cand_num)];
        else
            % Specrta available: Get candiate scores
            test_set = find (ismember (dt_inchi_mf_fp.mol_id, pathway_cmps_spec_id_best));
            assert (n_spec_best == numel (test_set));
            
            mf_test_set = dt_inchi_mf_fp.molecular_formula (test_set);
            assert (numel (unique (mf_test_set)) == 1)
            
            Y_C_test = CandidateSets (DataHandle (cand.cand), ...
                get_mf_corres (mf_test_set, cand.cand), 'ALL', fp_mask);
            
            fprintf ('\tCalculate kernels and scores ... ');
            tic;
            KX_list_train_test = cellfun(@(x) x(train_set, test_set), KX_list, 'UniformOutput', false);
            KX_list_test       = cellfun(@(x) x(test_set,  test_set), KX_list, 'UniformOutput', false);
            
            [scores, ~, hh] = Test_IOKR (KX_list_train_test, KX_list_test, ...
                train_model, Y_train, Y_C_test, param.iokr_param.center);
            fprintf ('%.3fs\n', toc);  
            
            fprintf ('\tWrite out scores for each spectra ... \n');
            tic;
            for idx_spec = 1:n_spec_best
                score = (scores{idx_spec} / sqrt(hh(idx_spec)))';
                id1 = Y_C_test.getCandidateSet (idx_spec, false, 'id')';
                
                [score_sorted, ix] = sort (score, 'descend');
                id1_sorted = id1(ix);
                
                score_list = table (id1_sorted, score_sorted);
                score_list.Properties.VariableNames = {'id1', 'score'};
                
                assert (any (Y_C_test.findExampleInCandidateSet (idx_spec, ...
                     pathway_cmps.inchikey_1(idx_mol))));       
                
                score_list_fn = sprintf ('%s/scoring_list=%s_spec=%s_id=%s.csv', ...
                    output_dir, ...
                    dt_inchi_mf_fp.molecular_formula{test_set(idx_spec)}, ...
                    dt_inchi_mf_fp.mol_id{test_set(idx_spec)}, ...
                    pathway_cmps.inchikey{idx_mol});
                fprintf ('\t\tScoring-file (%d/%d): %s\n', idx_spec, n_spec_best, ...
                    score_list_fn);
        
                writetable (score_list, score_list_fn, 'Delimiter', ',', ...
                    'WriteRowNames', false, 'WriteVariableNames', true')
                
                rank        = getRanksBasedOnScores(Y_C_test.getSubset(idx_spec), ...
                    pathway_cmps.inchikey_1(idx_mol), {score});
                cand_num    = Y_C_test.getCandidateSet(idx_spec, false, 'num');
                id          = pathway_cmps.inchikey(idx_mol);
                spec_id     = dt_inchi_mf_fp.mol_id(test_set(idx_spec));
                ranks_table = [ranks_table; table(id, rank, spec_id, cand_num)];
            end % for
            fprintf ('%.3fs\n', toc);      
            
            
            
%             ranks = [ranks; getRanksBasedOnScores(Y_C_test, ...
%                 repmat ({pathway_cmps.inchikey_1(idx_mol)}, [n_spec_best, 1]), scores)];
%             maxCandNum = max (maxCandNum, max (arrayfun (@(idx) Y_C_test.getCandidateSet (1, false, 'num'), ...
%                 1:Y_C_test.getNumberOfExamples())));
        end % if
    end % for
    
    
    rankPerc = getRankPerc (ranks_table.rank, max (ranks_table.cand_num));
    
    save ([output_dir, '/rank_perc.txt'], 'rankPerc', '-ascii');
    
    writetable ([ranks_table; ranks_table_no_spectra], [output_dir, '/ranks.csv'], 'Delimiter', ',', ...
            'WriteRowNames', false, 'WriteVariableNames', true');
    
    disp (rankPerc([1, 5, 10, 20]));
    

    

end % function
% 
%     mf_test_set = dt_inchi_mf_fp.molecular_formula (test_set);
%     Y_C_test = CandidateSets (DataHandle (cand.cand), get_mf_corres (mf_test_set, cand.cand), 'ALL', fp_mask);
%     
%     % Calulate the scores
%     fprintf ('Calculate scores ... ');
%     tic;
%     KX_list_train_test = cellfun(@(x) x(train_set, test_set), KX_list, 'UniformOutput', false);
%     KX_list_test       = cellfun(@(x) x(test_set,  test_set), KX_list, 'UniformOutput', false);
%     [scores, ~, hh] = Test_IOKR (KX_list_train_test, KX_list_test, ...
%             train_model, Y_train, Y_C_test, param.iokr_param.center);
%     fprintf ('%.3fs\n', toc);  
%         
%     % Create list of scores and save the result
%     %   scores_list=ML_spec=SPECTRANAME_id_=INCHIKEY.csv
%     spec_name_test_set = dt_inchi_mf_fp.mol_id(test_set);
%     id1_pathway_test_set = pathway_cmps.inchikey_1(locB_best(locB_best > 0));
%     id_pathway_test_set  = pathway_cmps.inchikey(locB_best(locB_best > 0));
%     ml_pathway_test_set  = pathway_cmps.ml_neutral(locB_best(locB_best > 0));
%     
    
%     fprintf ('Write out scores ...\n')
%     tic;
%     for idx = 1:Y_C_test.getNumberOfExamples()
%         fprintf ('\t%d/%d: ', idx, Y_C_test.getNumberOfExamples());
%         
%          score = (scores{idx} / sqrt (hh(idx)))';
%         
%         id1 = Y_C_test.getCandidateSet (idx, false, 'id')';
%         score_list = table (id1, score);
%         
%         assert (any (Y_C_test.findExampleInCandidateSet (idx, id1_pathway_test_set{idx})));
%         
%         score_list_fn = sprintf ('%s/scoring_list=%s_spec=%s_id=%s.csv', ...
%             output_dir, ml_pathway_test_set{idx}, spec_name_test_set{idx}, id_pathway_test_set{idx});
%         fprintf ('%s\n', score_list_fn);
%         
%         writetable (score_list, score_list_fn, 'Delimiter', ',', ...
%             'WriteRowNames', false, 'WriteVariableNames', true')
%     end % for
%     fprintf ('%.3fs\n', toc);
%     
%     fprintf ('Write out pseudo-scores ...\n');
%     has_no_spectra = true (size (pathway_cmps, 1), 1);
%     has_no_spectra(locB_best(locB_best > 0)) = false;
%     id_pathway_pseudo_set  = pathway_cmps.inchikey(has_no_spectra);
%     ml_pathway_pseudo_set  = pathway_cmps.ml_neutral(has_no_spectra);
%     
%     for idx = 1:sum (has_no_spectra)
%         score_list_fn = sprintf ('%s/scoring_list=%s_spec=NOSPECTRA_id=%s.csv', ...
%             output_dir, ml_pathway_pseudo_set{idx}, id_pathway_pseudo_set{idx});
%         fprintf ('%s\n', score_list_fn);
%         writetable (1, score_list_fn, 'Delimiter', ',', ...
%             'WriteRowNames', false, 'WriteVariableNames', true')
%     end % for
%     
%     % Calculate the top-k performance
%     ranks = getRanksBasedOnScores (Y_C_test, id1_pathway_test_set, scores);
%     maxCandNum = arrayfun (@(idx) Y_C_test.getCandidateSet (idx, 0, 'num'), 1:Y_C_test.getNumberOfExamples());
%     rankPerc = getRankPerc (ranks, maxCandNum);
    
%     pathway_cmps_spec_id_all = cellfun (@(c) strsplit (c, ';'), pathway_cmps.spec_id_neg_all, ...
%         'UniformOutput', false);
%     pathway_cmps_spec_id_all = 
%     
%     
%     
%     pathway_cmps_spec_id_best = cellfun (@(c) strsplit (c, ';'), pathway_cmps.spec_id_neg, ...
%         'UniformOutput', false);
%     
% 
%     
%     
%     
%     n_spec = size (dt_inchi_mf_fp, 1);
%     locB_all = zeros (n_spec, 1);
%     for idx = 1:n_spec
%         spec_id = dt_inchi_mf_fp.mol_id(idx);
%         tmp = cellfun (@(c) ismember (spec_id, c), pathway_cmps_spec_id_all);
%         
%         disp (tmp')
%         assert (sum (tmp) <= 1);
%         
%         if (any (tmp))
%             locB_all(idx) = find (tmp);
%         end % if
%     end % for
% %     [~, locB_all] = ismember (dt_inchi_mf_fp.mol_id, pathway_cmps_spec_id_all);
%     train_set = find (locB_all == 0); % All compounds not in the pathway
%     
%     % During the testing / scoring phase we use the spectra which has the
%     % best matching: 3D > 2D > 2D uncharged
%     locB_best = NaN (n_spec, 1);
%     for idx = 1:n_spec
%         spec_id = dt_inchi_mf_fp.mol_id(idx);
%         tmp = cellfun (@(c) ismember (spec_id, c), pathway_cmps_spec_id_best);
%         
%         assert (sum (tmp) <= 1);
%         
%         if (any (tmp))
%             locB_best(idx) = find (tmp);
%         end % if
%     end % for
% %     [~, locB_best] = ismember (dt_inchi_mf_fp.mol_id, pathway_cmps_spec_id_best);
%     test_set  = find (locB_best > 0);
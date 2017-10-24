function run_IOKR_for_candidate_scoring_special_course (input_dir, output_dir)
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
    
    % Load the pathway compounds:
    %   These compounds needs to be excluded from the training.
    impact_cmps = readtable ([input_dir, '/impact_cmps_cand.csv']);
    impact_cmps_mol_id_all = strsplit (strjoin(impact_cmps.mol_id_all, ';'), ';');
    impact_cmps_mol_id_all = impact_cmps_mol_id_all(cellfun(@(c) numel(c), impact_cmps_mol_id_all) > 0);
    
    % All the spectra which are associated with the any (3D, 2D, unchg 2D)
    % structures of the pathway compounds are excluded from the training.
    train_set = find (~ ismember (dt_inchi_mf_fp.mol_id, impact_cmps_mol_id_all));
   
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
    tic;
        load ([input_dir, '/cand_prep.mat']);
    fprintf ('%.3fs\n', toc);
    
    % Get the corresponding candidates set each spectra in the test set:
    %   For some pathway compounds we do have several spectra. We predict
    %   the scores for each candidate for each spectra, so that we can
    %   later for example average the scores belonging to different spectra
    %   measured from the same compound.
    eval = load ([input_dir, '/ind_eval.txt'], '-ascii');
    n_mol_pathway = size (impact_cmps, 1);
%     dt_inchi_mf_fp_eval = dt_inchi_mf_fp(eval, :);
    
    ranks_table = table;
    
    for idx_mol = 1:n_mol_pathway
%     for idx_mol = 170:n_mol_pathway
        fprintf ('Scoring: %d/%d\n', idx_mol, n_mol_pathway);
        
        impact_cmps_mol_id_best = strsplit (impact_cmps.mol_id_all{idx_mol}, ';');
        n_spec_best = numel (impact_cmps_mol_id_best);
        
        if ~ strcmp (impact_cmps_mol_id_best{1}, '')
            % Specrta available: Get candiate scores
            test_set = find (ismember (dt_inchi_mf_fp.mol_id, impact_cmps_mol_id_best));
            assert (n_spec_best == numel (test_set));
            
            % Restrict the score prediction to the evaluation set
            test_set = intersect (test_set, eval, 'stable');
            
            % Select randomly one spectra for the prediction
            if length (test_set) > 1
                test_set = randsample (test_set, 1);
            end % if
            n_spec_best = length (test_set);
            
            mf_test_set = dt_inchi_mf_fp.molecular_formula (test_set);
            assert (numel (unique (mf_test_set)) == 1)
            
            Y_C_test = CandidateSets (DataHandle (cand), ...
                get_mf_corres (mf_test_set, cand));
            
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
                     dt_inchi_mf_fp.inchi{test_set(idx_spec)})));       
                
                score_list_fn = sprintf ('%s/scoring_list=%s_spec=%s_id=%s.csv', ...
                    output_dir, ...
                    dt_inchi_mf_fp.molecular_formula{test_set(idx_spec)}, ...
                    dt_inchi_mf_fp.mol_id{test_set(idx_spec)}, ...
                    impact_cmps.inchikey{idx_mol});
                fprintf ('\t\tScoring-file (%d/%d): %s\n', idx_spec, n_spec_best, ...
                    score_list_fn);
        
                writetable (score_list, score_list_fn, 'Delimiter', ',', ...
                    'WriteRowNames', false, 'WriteVariableNames', true')
                
                rank        = getRanksBasedOnScores(Y_C_test.getSubset(idx_spec), ...
                    dt_inchi_mf_fp.inchi(test_set(idx_spec)), {score});
                cand_num    = Y_C_test.getCandidateSet(idx_spec, false, 'num');
                id          = impact_cmps.inchikey(idx_mol);
                mol_id      = dt_inchi_mf_fp.mol_id(test_set(idx_spec));
                ranks_table = [ranks_table; table(id, rank, mol_id, cand_num)];
            end % for
            fprintf ('%.3fs\n', toc);      
        end % if
    end % for
    
    
    rankPerc = getRankPerc (ranks_table.rank, max (ranks_table.cand_num));
    
    save ([output_dir, '/rank_perc.txt'], 'rankPerc', '-ascii');
    
    writetable (ranks_table, [output_dir, '/ranks.csv'], 'Delimiter', ',', ...
            'WriteRowNames', false, 'WriteVariableNames', true');
    
    disp (rankPerc([1, 5, 10, 20]));
   

end % function
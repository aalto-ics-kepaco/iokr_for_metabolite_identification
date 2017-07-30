%--------------------------------------------------------------
% Load and prepare data
%--------------------------------------------------------------

% inchi keys, molecular formulas, fingerprints
load ([inputDir '/compound_info.mat'], 'dt_inchi_mf_fp');
inchi = dt_inchi_mf_fp.inchi_key_1; 

% Extract fingerprints
Y = full (dt_inchi_mf_fp.fp_masked)';
[~,n] = size(Y);
mean_Y = mean (Y, 2);

% Candidates description
mf_corres = get_mf_corres (dt_inchi_mf_fp.molecular_formula, cand_transposed__);
Y_C       = CandidateSets (DataHandle (cand_transposed__), mf_corres);
assert (Y_C.getNumberOfExamples() == size (Y, 2))

%% ALL
[mean_all, cov_all] = Compute_cov_mean_feat (Y_C, mean_Y, true, true);

%% Perc
mean_perc = cell(0);
cov_perc  = cell(0);

ii = 1;
for perc = [1, 5, 10, 25, 50, 100];
    % Candidate selection
    param.data_param.selection_param = struct ( ...
            'strategy', 'random', 'perc', perc, 'inclExpCand', false);
    selec                            = getCandidateSelection (Y_C, inchi, ...
        param.data_param.selection_param);      
    Y_C.setSelectionsOfCandidateSets (selec);

    [mean_perc{ii}, cov_perc{ii}] = Compute_cov_mean_feat (Y_C, mean_Y, true, true);
    
    ii = ii + 1;
end % for 
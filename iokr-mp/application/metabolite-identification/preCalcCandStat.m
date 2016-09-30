function preCalcCandStat (cand)
    %% Load data
    tic;
    
    inputDir = '/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/';

    % Fingerprints
    Y = load (strcat (inputDir, '/fingerprints/fp.mat'));
    Y = full (Y.Y);

    % Identifier (inchis) for the training examples
    inchis = readtext (strcat (inputDir, '/inchi.txt'));
    
    % Pre-defined cross-validation folds to prevent to train using
    % molecular structures which are contained in the test
    ind_fold = load (strcat (inputDir, '/cv_ind.txt'));

    % Candidate sets 
    mf_corres = load (strcat (inputDir, '/candidates/matching_mf_train.txt'));
    if (nargin < 1)
        % For debuggin purposes we might want to pass in the candidate set.
        cand = load (strcat (inputDir, '/candidates/GNPS_cand_as_struct_transp.mat'));
        cand = cand.cand_struct;
    end % if
    Y_C = CandidateSets (DataHandle (cand), mf_corres);

    clear cand;

    fprintf ('Loading data took %.3fs\n', toc);

    %% Set parameter
    param = struct ();
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (param, ...
        {'debug_param', 'opt_param', 'mp_iokr_param', 'data_param'});

    % Cross-validation settings
    cv_param = struct ('outer', struct ('type', 'fixed', 'cvInd', ind_fold), ...
                       'inner', struct ('nFolds', param.opt_param.nInnerFolds));                       
    param.opt_param.cv_param = cv_param;

    % Candidate selection settings
    param.data_param.selection_param = struct ('strategy', 'all', 'inclExpCand', true);

    %% Pre-calculate statistics
    tic;
    getPreCalcCandStat_feat (Y, Y_C, inchis, param, inputDir);
    fprintf ('Loading / pre-calculating of the candidate statistics took %.3fs\n', toc);
end % if
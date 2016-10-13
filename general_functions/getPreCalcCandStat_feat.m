function matObj = getPreCalcCandStat_feat (Y, Y_C, inchis, param, inOutDir, do_not_modify_Y_C)
%% PRECALCCANDSTAD_FEAT preculates the candidate statistics given a set of parameters
%
%       PARAM.TMP_PARAM.CV must correspont to GETCVINDICES (PARAM.DATA_PARAM.CV)
%       PARAM.TMP_PARAM.SELEC must correspont to GETCANDIDATESELECTION (Y_C, inchis, PARAM.TMP_PARAM.SELEC)
%
%    INPUTS:
%       Y_C     CandidateSets object. The selection must be set
%               corresponding to PARAM.DATA_PARAM.SELECTION_PARAM.
%
%    TODO:
%       INCHIS should not be a parameter of this function. However, in order
%       to prevent miss use (by me ;)) it will be passed into the function
%       for the moment. FIXME: INCHI should go as identifier into Y_C!
%
%       The parameter used for pre-calculation must be associated with the
%       stored statistics. Those statistics are stored using a filename
%       which is just a hash-value based on the underlying parameters. We
%       need a file, database, etc. which associates a hash with a set of
%       paramters.
    if (nargin < 6)
        warning ('Please fix this in the future: do_not_modify_Y_C');
        do_not_modify_Y_C = false;
    end % if

    sw_stats    = StopWatch ('Statistics outer fold (WALL)');
    sw_stats_cv = StopWatch ('Statistics inner fold (WALL)');
    sw_fold     = StopWatch ('One outer fold (WALL)');
    
    swc_stats    = StopWatchCPUTime ('Statistics outer fold (CPU)');
    swc_stats_cv = StopWatchCPUTime ('Statistics inner fold (CPU)');
    swc_fold     = StopWatchCPUTime ('One outer fold (CPU)');

    % If precaclulated statistics are used, the following things have
    % to be loaded:
    %   1.1) cross-validation data splitting
    %   1.2) candidate selection 
    %   1.3) repetition, used for the random selection
    %   1.4) centering setting
    %   1.5) Run in debug mode? (if yes: use settings of the debug mode)
    %   --> covariance matrices & mean vectors
    % Just to be sure: Everything will be in ONE file. 

    % Use hash to find file matching the parameters
    param_struct = struct (                                  ...
        'cv_param',        param.data_param.cv_param,        ...
        'selection_param', param.data_param.selection_param, ...
        'repetition',      param.data_param.repetition,      ...
        'center',          param.mp_iokr_param.center,       ...
        'is_debug_mode',   param.debug_param.isDebugMode);
    
    if (param_struct.is_debug_mode)
        param_struct.('random_seed') = param.debug_param.randomSeed;
        param_struct.('n_debug_set') = param.debug_param.n_debug_set;
    else
        param_struct.('random_seed') = NaN;
        param_struct.('n_debug_set') = false(0);
    end % if
    
    % Create filename based on the settings
    statHash = DataHash (param_struct);
    statFn = strcat (inOutDir, '/', statHash, '.mat');
    
    % Load respectively pre-calculated the statistics
    if (~ exist (statFn, 'file'))
        warning ('No pre-calculated statistics available. %s: No such file.', ...
            statFn);
        disp ('Pre-calculate statistics');

        try 
            cv    = getCVIndices (param.data_param.cv_param);
            selec = getCandidateSelection (Y_C, inchis, param.data_param.selection_param);
            
            Y_C.setSelectionsOfCandidateSets (selec);
            
            % Create a new statistic file
            matObj = matfile (statFn, 'Writable', true);

            matObj.cv         = cv;
            matObj.selec      = selec;
            
            matObj.repetition = param.data_param.repetition;
            matObj.center     = param.mp_iokr_param.center;

            % Temporay variables
            ker_center  = param.mp_iokr_param.center;
            nOuterFolds = cv.outer.NumTestSets;
            nInnerFolds = cv.inner{1}.NumTestSets;
            
            matObj.stats    = repmat (struct ('Mean_Psi_C_train', [], 'Cov_Psi_C_train', []), ...
                [1, nOuterFolds]);
            matObj.stats_cv = repmat (struct ('Mean_Psi_C_train_cv', [], 'Cov_Psi_C_train_cv', [], ...
                'Mean_Psi_C_test_cv', [], 'Cov_Psi_C_test_cv', []), [nInnerFolds, nOuterFolds]);
            

            assert (all (cellfun (@(c) c.NumTestSets, cv.inner) == nInnerFolds), ...
                'The nunber of inner folds must be equal for all outer folds.');
            
            for iOutFold = 1:nOuterFolds
                sw_fold.start();
                swc_fold.start();
                fprintf ('Outer fold: %d/%d\n', iOutFold, nOuterFolds);
                
                train_set = training_my (cv.outer, iOutFold);
                                
                Y_train = Y(:, train_set);
                
                mean_Y_train = mean (Y_train, 2);
                
                Y_C_train = Y_C.getSubset (train_set);

                sw_stats.start();
                swc_stats.start();
                % Train statistics
                [Mean_Psi_C_train, Cov_Psi_C_train] = ...
                    Compute_cov_mean_feat (Y_C_train, mean_Y_train, ker_center, true);   
                sw_stats.stop();
                swc_stats.stop();
                sw_stats.showAvgTime();
                swc_stats.showAvgTime();
                
                matObj.stats(1, iOutFold) = struct (      ...
                    'Mean_Psi_C_train', Mean_Psi_C_train, ...
                    'Cov_Psi_C_train', Cov_Psi_C_train);

                for iInFold = 1:nInnerFolds
                    fprintf ('Inner fold: %d/%d\n', iInFold, nInnerFolds);
                    
                    train_set_cv = training_my (cv.inner{iOutFold}, iInFold);
                    test_set_cv  = test_my (cv.inner{iOutFold}, iInFold);

                    Y_train_cv = Y_train(:, train_set_cv);
                    Y_test_cv  = Y_train(:, test_set_cv);

                    mean_Y_train_cv = mean (Y_train_cv, 2);
                    mean_Y_test_cv  = mean (Y_test_cv, 2);

                    Y_C_train_cv = Y_C_train.getSubset (train_set_cv);
                    Y_C_test_cv  = Y_C_train.getSubset (test_set_cv);

                    sw_stats_cv.start();
                    swc_stats_cv.start();
                    % Train / Test (cv) statistics
                    [Mean_Psi_C_train_cv, Cov_Psi_C_train_cv] = Compute_cov_mean_feat (Y_C_train_cv, mean_Y_train_cv, ker_center, true);
                    [Mean_Psi_C_test_cv, Cov_Psi_C_test_cv]   = Compute_cov_mean_feat (Y_C_test_cv, mean_Y_test_cv, ker_center, true);
                    sw_stats_cv.stop();
                    swc_stats_cv.stop();
                    sw_stats_cv.showAvgTime();
                    swc_stats_cv.showAvgTime();

                    if ((iOutFold == 1) && (iInFold == 1))
                    matObj.stats_cv = struct (   ...
                            'Mean_Psi_C_train_cv', Mean_Psi_C_train_cv, ...
                            'Mean_Psi_C_test_cv',  Mean_Psi_C_test_cv,  ...
                            'Cov_Psi_C_train_cv',  Cov_Psi_C_train_cv,  ...
                            'Cov_Psi_C_test_cv',   Cov_Psi_C_test_cv);    
                    else
                        matObj.stats_cv(iInFold, iOutFold) = struct (   ...
                            'Mean_Psi_C_train_cv', Mean_Psi_C_train_cv, ...
                            'Mean_Psi_C_test_cv',  Mean_Psi_C_test_cv,  ...
                            'Cov_Psi_C_train_cv',  Cov_Psi_C_train_cv,  ...
                            'Cov_Psi_C_test_cv',   Cov_Psi_C_test_cv);                         
                    end % if
                end % for
                
                sw_fold.stop();
                swc_fold.stop();
                sw_fold.showAvgTime();
                swc_fold.showAvgTime();
                
                fprintf ('Remaining time (estimate): %.3fh\n', ...
                    ((nOuterFolds - iOutFold) * sw_fold.getAverageTime()) / 3600);
            end % for
                        
            clear matObj;
        catch me
            % If anything goes wrong:
            % ... the statistic file should be removed
            delete (statFn);
            % ... the candidate selection should be reset
            Y_C.resetSelectionToAllCandidates();
            % This leaves the world as it was. 

            rethrow (me);
        end % try      
    end % if
    
    if (nargout > 0)
        disp (['Load pre-caculated statistics from file: ', statFn]);
        matObj = matfile (statFn);
        
        if (do_not_modify_Y_C) ; return ; end % if
        
        Y_C.setSelectionsOfCandidateSets (matObj.selec);
    end % if
end % function
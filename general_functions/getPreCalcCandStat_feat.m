function matObj = getPreCalcCandStat_feat (Y, Y_C, inchis, param, inOutDir)
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
%       for the moment. Btw: INCHI should go as identifier into Y_C!
%
%       The parameter used for pre-calculation must be associated with the
%       stored statistics. Those statistics are stored using a filename
%       which is just a hash-value based on the underlying parameters. We
%       need a file, database, etc. which associates a hash with a set of
%       paramters.

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
        'cv_param',        param.opt_param.cv_param,         ...
        'selection_param', param.data_param.selection_param, ...
        'repetition',      param.data_param.repetition,      ...
        'center',          param.mp_iokr_param.center,       ...
        'is_debug_mode',   param.debug_param.isDebugMode);
    
    if (param_struct.is_debug_mode)
        param_struct.('random_seed') = param.debug_param.isDebugMode;
        param_struct.('n_debug_set') = param.debug_param.n_debug_set;
    else
        param_struct.('random_seed') = NaN;
        param_struct.('n_debug_set') = false(0);
    end % if
    
    % Create filename based on the settings
    statHash = DataHash (param_struct);
    statFn = strcat (inOutDir, '/pre_calculated_stats/', statHash, '.mat');
    
    cv    = getCVIndices (param.opt_param.cv_param);
    selec = getCandidateSelection (Y_C, inchis, param.data_param.selection_param);
    
    % Load respectively pre-calculated the statistics
    if (~ exist (statFn, 'file'))
        warning ('No pre-calculated statistics available. %s: No such file.', ...
            statFn);
        disp ('Pre-calculate statistics');

        try 
            Y_C.setSelectionsOfCandidateSets (selec);
            
            % Create a new statistic file
            matObj = matfile (statFn, 'Writable', true);

            matObj.cv         = cv;
            matObj.selec      = selec;
            
            matObj.repetition = param.data_param.repetition;
            matObj.center     = param.mp_iokr_param.center;
            % Dimension: [1 x param.data_param.cv.outer.NumTestSets]
            matObj.stats      = struct ([]);
            % Dimension: [param.data_param.cv.inner.NumTestSets x cv.outer.NumTestSets]
            matObj.stats_cv   = struct ([]);

            for i = 1:cv.outer.NumTestSets
                fprintf ('Outer fold: %d/%d\n', i, cv.outer.NumTestSets);
                
                train_set = training_my (cv.outer, i);
                                
                Y_train = Y(:, train_set);
                
                mean_Y_train = mean (Y_train, 2);
                
                Y_C_train = Y_C.getSubset (train_set);

                tic;
                [Mean_Psi_C_train, Cov_Psi_C_train] = Compute_cov_mean_feat ( ...
                    Y_C_train, mean_Y_train, param.mp_iokr_param.center);     
                toc;

                if (i == 1)
                    matObj.stats = struct ('Mean_Psi_C_train', Mean_Psi_C_train, ...
                        'Cov_Psi_C_train', Cov_Psi_C_train);
                else
                    matObj.stats(1, i) = struct ('Mean_Psi_C_train', Mean_Psi_C_train, ...
                        'Cov_Psi_C_train', Cov_Psi_C_train);
                end % if
                
                clear Mean_Psi_C_train Cov_Psi_C_train;

                for j = 1:cv.inner{i}.NumTestSets
                    fprintf ('Inner fold: %d/%d\n', j, cv.inner{i}.NumTestSets);
                    
                    train_set_cv = training_my (cv.inner{i}, j);
                    test_set_cv  = test_my (cv.inner{i}, j);

                    Y_train_cv = Y_train(:, train_set_cv);
                    Y_test_cv  = Y_train(:, test_set_cv);

                    mean_Y_train_cv = mean (Y_train_cv, 2);
                    mean_Y_test_cv  = mean (Y_test_cv, 2);

                    Y_C_train_cv = Y_C_train.getSubset (train_set_cv);
                    Y_C_test_cv  = Y_C_train.getSubset (test_set_cv);

                    tic;
                    [Mean_Psi_C_train_cv, Cov_Psi_C_train_cv] = Compute_cov_mean_feat ( ...
                        Y_C_train_cv, mean_Y_train_cv, param.mp_iokr_param.center);     
                    toc;
                    tic;
                    [Mean_Psi_C_test_cv, Cov_Psi_C_test_cv] = Compute_cov_mean_feat ( ...
                        Y_C_test_cv, mean_Y_test_cv, param.mp_iokr_param.center);
                    toc;

                    if ((i == 1) && (j == 1))
                        matObj.stats_cv = struct ('Mean_Psi_C_train_cv', Mean_Psi_C_train_cv, ...
                                                  'Mean_Psi_C_test_cv',  Mean_Psi_C_test_cv,  ...
                                                  'Cov_Psi_C_train_cv',  Cov_Psi_C_train_cv,  ...
                                                  'Cov_Psi_C_test_cv',   Cov_Psi_C_test_cv);
                    else
                        matObj.stats_cv(j, i) = struct ('Mean_Psi_C_train_cv', Mean_Psi_C_train_cv, ...
                                                        'Mean_Psi_C_test_cv',  Mean_Psi_C_test_cv,  ...
                                                        'Cov_Psi_C_train_cv',  Cov_Psi_C_train_cv,  ...
                                                        'Cov_Psi_C_test_cv',   Cov_Psi_C_test_cv);                                                       
                    end % if

                    clear Mean_Psi_C_train_cv Cov_Psi_C_train_cv ...
                        train_set_cv test_set_cv ...
                        Y_train_cv mean_Y_train_cv Y_C_train_cv ...
                        Y_test_cv mean_Y_test_cv Y_C_test_cv;
                end % for
            end % for
            
            % Check post-condition
            assert (all (size (matObj.stats) == [1, cv.outer.NumTestSets]));
            cellfun (@(c) assert ( ...
                all (size (matObj.stats_cv) == [c.NumTestSets, cv.outer.NumTestSets])), cv.inner);
                        
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
    
    disp (['Load pre-caculated statistics from file: ', statFn]);

    matObj = matfile (statFn);
    
    Y_C.setSelectionsOfCandidateSets (matObj.selec);
end % function
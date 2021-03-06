function cv = getCVIndices (cvParam, modeStr, directory, overwrite)
%% GETCVINDICES Function which pre-computes the train/test indices for CV
%    cv = GETCVINDICES (CVPARAM) creates a partition according to the 
%    cross-validation setting CVPARAM.
%
%    cv = GETCVINDICES (CVPARAM, 'cache', DIRECTORY) creates a partition
%    according to the provided parameter structure and stores this
%    partition under in DIRECTORY. If the file already exists a warning is 
%    provided and nothing is written out. 
%
%    cv = GETCVINDICES (CVPARAM, 'cache', DIRECTORY, OVERWRITE) if the
%    binary value OVERWRITE is set to true, than a possibly existing
%    partition will be overwritten.
%
%    cv = GETCVINDICES (CVPARAM, 'load', DIRECTORY) loads a partition from 
%    the DIRECTORY.
% 
%    INPUTS
%       cvParam             Parameter structure to build up the nested-cv.
%    * Example 1.1:
%       cvParam = struct ('nObservations', nObservations,    ...
%                   'outer', struct ('type', 'random',       ...
%                                    'nFolds', nOuterFolds), ... 
%                   'inner', struct ('nFolds', nInnerFolds))
%
%       Outer and inner folds are chosen randomly. The number of outer and
%       inner folds must be provided as well as the number of examples 
%       NOBSERVATIONS.
%
%    * Example 1.2:
%       cvParam = struct ('nObservations', nObservations, ...
%                   'outer', struct ('type', 'random',    ...
%                                    'nFolds', nOuterFolds))
% 
%       Same as example 1.1, but no inner cross-validation folds will be
%       created.
%
%    * Example 1.3:
%       cvParam = struct ('nObservations', nObservations, ...
%                   'randomSeed', fixedRandomSeed,        ...
%                   'outer', struct ('type', 'random',    ...
%                                    'nFolds', nOuterFolds))
%
%       'fixedRandomSeed' is an integer value which is used as random seed
%       for the pseudo random generator (PRG). The current state of the
%       PRG is stored and restored after the outer and inner cv-folds have
%       been chosen randomly using the fixed seed.
%
%    * Example 2.1:
%       cvParam = struct (                            ...
%                   'outer', struct ('type', 'fixed', ... 
%                                    'cvInd', cvInd), ... 
%                   'inner', struct ('nFolds', nInnerFolds))
%
%       Outer folds are fixed using a vector of fold indices. The number of
%       outer folds is determined using 'numel (unique (cvInd))' and the
%       number of training examples using 'numel (cvInd)'. Each element of
%       'cvInd(i)' is treat as the cv-fold of training example i. The inner
%       folds are chosen randomly and the number of inner folds should be
%       provided.
%
%    * Example 2.2:
%       cvParam = struct (                            ...
%                   'outer', struct ('type', 'fixed', ... 
%                                    'cvInd', cvInd))
%
%       directory           Output/input directory to store/load the
%                           candidate selection.
%       overwrite           Binary indicating whether an exsiting partition
%                           should be overwritten.
%
%       Same as example 2.1, but no inner cross-validation folds will be
%       created.
%
%    OUTPUTS:
%       cv    Struct with the following structure:
%             * outer: Depending on PARAM.OUTER.TYPE: 
%               'random' --> Object of class cvpartition
%               'fixed'  --> Struct:
%                            * training:    (nTrain x nOuterFolds) logical
%                                           matrix. (:, i) corresponds to
%                                           the selection of the training
%                                           examples for outer fold i.
%                            * test:        (nTrain x nOuterFolds) logical
%                                           matrix. (:, i) corresponds to
%                                           the selection of the test
%                                           examples for outer fold i.
%             If requested by provided PARAM.INNER structure. 
%             * inner: (nOuterFolds x 1) cell array containing objects of 
%                      class cvpartition.
%
%    NOTE: The class cvpartitions cannot work with pre-defined cv-folds.
%          However, to provide the following functionalities:
%
%          foldSelec = training_my (cv.outer, i)
%          foldSelec = test_my (cv.outer, i)
%
%          foldSelec = training_my (cv.inner{i}, j)
%          foldSelec = test_my (cv.inner{i}, j)
%           
%          The functions TRAINING_MY and TEST_MY has been implemented.
%          Those function can work with object of class cvpartition and
%          simple structures created by this function.
%
%    See also CVPARTITION, TRAINING_MY and TEST_MY.
    if (nargin < 1)
        error ('getCVIndices:InvalidArgument', ...
            'Not enough arguments.');
    end % if
    
    [isValid, errStr] = validateParamStructure (cvParam);
    if (~ isValid) 
        error ('getCVIndices:InvalidArgument', 'Message: %s', errStr);
    end % if
    
    if (nargin < 2)
        modeStr = 'normal';
    else % nargin >= 2
        if (nargin < 3)
            error ('getCVIndices:InvalidArgument', ...
                'Working mode is "%s", but not directory path is specified.', modeStr);
        end % if
        [isValid, errStr] = validateDirectory (directory);
        if (~ isValid) 
            error ('getCVIndices:InvalidArgument', 'Message: %s', errStr);
        end % if
        
        if (nargin < 4)
            overwrite = false;
        end % if
    end % if
    
    % Calculate the cross-validation partition
    if (ismember (modeStr, {'normal', 'cache'}))
        switch (cvParam.outer.type)
            case 'random'
                nFoldsOuter = cvParam.outer.nFolds;                
                nObservations = cvParam.nObservations;
                
                % We can explicitly set the random seed to a desired, e.g.
                % fixed value, in order to ensure the cross-validation
                % folds to be the same across experiments.
                if (ismember ('randomSeed', fieldnames (cvParam)))
                    % Store the current random-state.
                    oldRandomState = rng;
                    % Set the current random-seed to a fixed one.
                    rng (cvParam.randomSeed);
                    
                    resetRandomSeed = true;
                else 
                    resetRandomSeed = false;
                end % if
                
                cv = struct ('outer', {});
                if (nFoldsOuter > 1)
                    cv(1).outer = cvpartition (nObservations, 'KFold', nFoldsOuter);
                else
                    cv(1).outer = cvpartition (nObservations, 'Resubstitution');
                end % if
                
                if (ismember ('inner', fieldnames (cvParam)))
                    nFoldsInner = cvParam.inner.nFolds;
                    
                    cv.inner = cell (nFoldsOuter, 1);
                    for i = 1:nFoldsOuter
                        cv(1).inner{i} = cvpartition (sum (cv.outer.training (i)), 'KFold', nFoldsInner);
                    end % for 
                else 
                    nFoldsInner = NaN;
                end % if
                
                if (resetRandomSeed)
                    % Restore random-state
                    rng (oldRandomState);
                end % if
            case 'fixed'
                nFoldsOuter = numel (unique (cvParam.outer.cvInd));              
                nObservations = numel (cvParam.outer.cvInd);
                
                cv = struct ( ...
                    'outer', struct ('training', false (nObservations, nFoldsOuter), ...
                                     'test', false (nObservations, nFoldsOuter),     ...
                                     'NumTestSets', nFoldsOuter));
                for i = 1:nFoldsOuter 
                    cv.outer.test(:, i)     = (cvParam.outer.cvInd == i);
                    cv.outer.training(:, i) = ~ cv.outer.test(:, i);
                end % for
                
                if (ismember ('inner', fieldnames (cvParam)))
                    nFoldsInner = cvParam.inner.nFolds; 
                    
                    cv.inner = cell(nFoldsOuter, 1); 
                    for i = 1:nFoldsOuter 
                        cv.inner{i} = cvpartition (sum (cv.outer.training (:, i)), 'KFold', nFoldsInner);
                    end % for
                else
                    nFoldsInner = NaN;
                end % if
            otherwise 
                error ('getCVIndices:InvalidArgument', ...
                    '"%s" is not a valid partition mode for the outer fold. Only "random" and "fixed" allowed.', ...
                    cvParam.outer.type);
        end % switch
        
        % Write out selection if needed
        if (strcmp (modeStr, 'cache'))
            cvFn = strcat (directory, '/cvPartition', statisticSettings2Str ( struct ( ...
                'foldOuter',    nFoldsOuter, ...
                'foldInner',    nFoldsInner, ...
                'fixedOuterCV', strcmp (cvParam.outer.type, 'fixed'))), '.mat');            
            
            if ((~ overwrite) && (exsit (cvFn, 'file')))
                warning ('%s already exist. Nothing will be written out.', ...
                    cvFn);
            else
                save (cvFn, 'cv', '-v7.3');
            end % if      
        end % if
    elseif (strcmp (modeStr, 'load'))
        switch (cvParam.outer.type)
            case 'random'
                nFoldsOuter = cvParam.outer.nFolds;
            case 'fixed'
                nFoldsOuter = numel (unique (cvParam.outer.cvInd));
            otherwise
                error ('getCVIndices:InvalidArgument', ...
                    '"%s" is not a valid partition mode for the outer fold. Only "random" and "fixed" allowed.', ...
                    cvParam.outer.type);
        end % if
        
        if (ismember ('inner', fieldnames (cvParam)))
            nFoldsInner = cvParam.inner.nFolds;      
        else
            nFoldsInner = NaN;
        end % if
        
        cvFn = strcat (directory, '/cvPartition', statisticSettings2Str ( struct ( ...
            'foldOuter',    nFoldsOuter, ...
            'foldInner',    nFoldsInner, ...
            'fixedOuterCV', strcmp (cvParam.outer.type, 'fixed'))), '.mat');
        
        tmp = load (cvFn);
        cv = tmp.cv;
    else
        error ('getCVIndices:InvalidArgument', ...
            '"%s" is not a valid working-mode. Leave parameter empty or set it to "load" or "cache".', ...
            modeStr);
    end % if   
end % function

% Function to validate the cross-validation parameter structure
function [isValid, errStr] = validateParamStructure (param)
    isValid = true;
    errStr = '';
    
    if (~ isa (param, 'struct'))
        isValid = false;
        errStr = 'Parameters must be provided as struct.';
        
        return;
    end % if
    if (~ ismember ('outer', fieldnames (param)))
        isValid = false;
        errStr = 'Parameter struct must contain the field "outer".';
        
        return;
    end % if
    
%     warning ('To Eric: Just write the input parameter validation ;).');
end % function

% Function to validate the working-directory
function [isValid, errStr] = validateDirectory (directory)
    isValid = true;
    errStr = '';

    if (~ ischar (directory))
        errStr = sprintf ('Directory path must be provided as string. Provided class = %s.', ...
            class (directory));
        return; 
    end % if
    if (~ exist (directory, 'dir'))
        isValid = false;
        errStr = sprintf ('"%s": No such directory.', directory);

        return;
    end % if
end % function
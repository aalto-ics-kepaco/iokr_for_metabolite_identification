function [selec, debugVec] = getCandidateSelection (candidateSets, id, param, ...
    modeStr, directory, overwrite)
%% GETCANDIDATESELECTION Function to select subset of candidates for MP
%    SELEC = GETCANDIDATESELECTION (CANDIDATESETS, ID, PARAM) returns the 
%    candidate selection. ID contains the identifier for each example, 
%    which is used to find the candidate corresponding to the example.
%    PARAM is structure defining the candidate selection strategy.
%     
%    SELEC = GETCANDIDATESELECTION (...,  'cache', DIRECTORY) The candidate 
%    selection is stored in DIRECTORY. The filename is determined using 
%    PARAM. If the file already exists a warning is provided and _nothing_ 
%    is written out. 
%
%    SELEC = GETCANDIDATESELECTION (..., 'cache', DIRECTORY, OVERWRITE) The 
%    binary OVERWRITE indicates, whether an exisiting candidate selection 
%    should be overwritten. 
% 
%    SELEC = GETCANDIDATESELECTION (CANDIDATESETS, [], PARAM, 'load', DIRECTORY) 
%    loads a candidate selection from the DIRECTORY. The filename is 
%    determined using PARAM. The CANDIDATESETS are updated to the loaded
%    selelection
%
%    INPUTS:
%       candidateSets           Object of class CANDIDATESETS. 
%       id                      (l x 1)-dimensional cell-array storing the
%                               identifier of each example
%       param                   Struct containing parameter for the 
%                               candidate selection.
%           * strategy:    String, either 'all' or 'random'
%           * inclExpCand: Binary indicating, whether the correct candidate 
%                          for a training example should be included in the 
%                          candidate selection or not. 
%           * perc:        Percentage of candidate _randomly_ chosen from the 
%                          candidate set for the selection (only for 'random')
%       directory               Output/input directory to store/load the
%                               candidate selection.
%       overwrite               Binary indicating whether an exsiting
%                               selection should be overwritten.
%
%    OUTPUTS:
%       selec                   (l x 1)-dimensional cell-array storing the
%                               logical-vectors representing the candidate
%                               selection. 
    if (nargin < 3)
        error ('getCandidateSelection:InvalidArgument', ...
            'Not enough arguments.');
    end % if
    
    if (~ isa (candidateSets, 'CandidateSets'))
        error ('getCandidateSelection:InvalidArgument', ...
            'Candidate sets must be provided using the CANDIDATESETS class.');
    end % if
    [isValid, errStr] = validateParamStructure (param);
    if (~ isValid) 
        error ('getCandidateSelection:InvalidArgument', 'Message: %s', errStr);
    end % if
    
    if (nargin < 4)
        modeStr = 'normal';
    else % nargin >= 4
        if (nargin < 5)
            error ('getCandidateSelection:InvalidArgument', ...
                'Working mode is "%s", but not directory path is specified.', modeStr);
        end % if
        [isValid, errStr] = validateDirectory (directory);
        if (~ isValid) 
            error ('getCandidateSelection:InvalidArgument', 'Message: %s', errStr);
        end % if
        
        if (nargin < 6)
            overwrite = false;
        end % if
    end % if
    
    if (nargout > 1)
        debugVec = true (numel (id), 1);
    end % if
    
    % Perform / load selection
    if (ismember (modeStr, {'normal', 'cache'}))
        [isValid, errStr] = validateIdStructure (id, candidateSets);
        if (~ isValid) 
            error ('getCandidateSelection:InvalidArgument', 'Message: %s', errStr);
        end % if
        
        % Calculate the candidate selection
        switch (param.strategy)
            case 'all'
                if (param.inclExpCand)
                    selec = candidateSets.createSelectionOfAllCandidates();
                else
                    numExamples = candidateSets.getNumberOfExamples();
                    selec = cell (numExamples, 1);

                    for i = 1:numExamples;
                        selec{i} = candidateSets.findExampleInCandidateSet (i, id(i));
                        
                        % We want to select everything _but_ the examples
                        % candidate
                        selec{i} = ~ selec{i};
                    end % for
                end % if
            case 'random'
                numExamples = candidateSets.getNumberOfExamples();
                selec = cell (numExamples, 1);
                
                % Save some statistics
                numSelecCand = 0;
                numCand      = 0;

                for i = 1:numExamples;
                    numCandForExamples = candidateSets.getCandidateSet (i, 0, 'num');
                    if (isnan (numCandForExamples))
                        % There is no candidate set
                        selec{i} = NaN;
                        
                        if (nargout > 1)
                            debugVec(i, :) = false;
                        end % if
                        
                        continue;
                    end % if
                    
                    exampleSelection = candidateSets.findExampleInCandidateSet (i, id(i));
                    if (~ any (exampleSelection))      
                        % The example is _not_ in the candidate set
                        r = rand (candidateSets.getCandidateSet (i, 0, 'num'), 1);
                        if (nargout > 1)
                            debugVec(i) = false;
                        end % if
                    else
                        % The example is in the candidate set
                        r = rand (candidateSets.getCandidateSet (i, 0, 'num') - 1, 1);
                    end %
                    % Add the selected candidates
                    selec{i}(~ exampleSelection) = r < (param.perc / 100);
                    % Add the example's candidate if needed
                    selec{i}(exampleSelection)   = logical (param.inclExpCand);
                    
                    numCand      = numCand + numCandForExamples;
                    numSelecCand = numSelecCand + sum (selec{i});
                end % for
                
                fprintf ('Number of all the candidates (NOT unique ones): %d\n', numCand);
                fprintf ('Number of selected candidates: %d\n', numSelecCand);
                fprintf ('=> %.3f%%\n', (100 * numSelecCand) / numCand);
            otherwise
                error ('getCandidateSelection:InvalidArgument', ...
                    'The parameter-structure should not contain any other strategy than "all" and "random"');   
        end % switch
        
        % Some assertions
        assert (all (cellfun (@(c) size(c, 1), selec) == 1));
        
        numCandPerSet = arrayfun (@(idx) candidateSets.getCandidateSet (idx, 0, 'num'), 1:candidateSets.getNumberOfExamples());
        numSelectionsPerSet = cellfun (@(c) size(c, 2), selec);
        assert (all (numSelectionsPerSet(~ isnan (numCandPerSet)) == numCandPerSet(~ isnan (numCandPerSet))'));
        
        % Write out selection if needed
        if (strcmp (modeStr, 'cache'))
            selecFn = strcat (directory, '/candidateSelection', statisticSettings2Str (param), '.mat');
            
            if ((~ overwrite) && (exsit (selecFn, 'file')))
                warning ('%s already exist. Nothing will be written out.', ...
                    selecFn);
            else
                save (selecFn, 'selec', '-v7.3');
            end % if      
        end % if
    elseif (strcmp (modeStr, 'load'))
        selecFn = strcat (directory, '/candidateSelection', statisticSettings2Str (param), '.mat');
        
        tmp = load (selecFn);
        selec = tmp.selec;
    else 
        error ('getCandidateSelection:InvalidArgument', ...
            '"%s" is not a valid working-mode. Leave parameter empty or set it to "load" or "cache".', ...
            modeStr);
    end % if   
end % function 

% Function to validate the id-structure
function [isValid, errStr] = validateIdStructure (id, candidateSets)
    isValid = true;
    errStr = '';
    
    if (all (size (id) > 1))
        isValid = false;
        errStr = 'Identifier must be provided in a vector.';
        
        return;
    end % if
    if (~ (numel (id) == candidateSets.getNumberOfExamples()))
        isValid = false;
        errStr = 'Number of identifiers must match the number of examples associated with the candidate sets.';
        
        return;
    end % if
end % function

% Function to validate the param-structure
function [isValid, errStr] = validateParamStructure (param)
    requieredFields = {'strategy', 'inclExpCand'};
    %validFields = struct ('strategy', {{'all', 'random'}});

    isValid = true;
    errStr = '';
    
    if (~ all (ismember (requieredFields, fieldnames (param))))
        isValid = false;
        errStr = sprintf ('Parameter structure must contain the following fields: %s\n', ...
            sprintf ('%s, ', requieredFields{:}));
        
        return;
    end % if
    
    switch (param.strategy)
        case 'random'
            if (~ (ismember ('perc', fieldnames (param))))
                isValid = false;
                errStr = 'If the selection strategy is random, than parameter structure must contain the field "perc"';
                
                return;
            end % if
    end % switch
    
%    for fn = requieredFields
%        switch fn{1}
%            case 'strategy'
%                if (~ ismember (param.(fn{1}), validFields.strategy))
%                    isValid = false;
%                    errStr = sprintf ('Only the following selection strategies are valid: %s\n', ...
%                        sprintf ('%s, ', validFields.strategy));
%                    
%                    return;
%                end % if
%            case 'inclExpCand'
%                if (~ islogical (param.(fn{1})))
%                    isValid = false;
%                    errStr = sprintf ('"param.inclExpCand" must be a logical. Provided class = %s.', ...
%                        class (param.(fn{1})));
%                    
%                    return;
%                end % if
%        end % switch
%    end % for
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
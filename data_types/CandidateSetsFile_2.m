classdef CandidateSetsFile_2 < handle
    properties(SetAccess = private, GetAccess = public)
        data_dir_;
        lut_;
        selec_feature_;
        n_features_;
    end % properties: write private, read public
    
    properties(Constant, Access = private)
        candidate_set_fieldnames = struct( ...
            'data', 'fingerprints', ...        
            'id', 'inchikey2D', ...
            'num', 'numberOfCandidates', ...
            'inchi2D', 'inchi2D');
    end % properties: constant, private
    
    methods(Access = public)
        function obj = CandidateSetsFile_2(data_dir, lut, selec_feature)
        %% CANDIDATESETSFILE Class-constructor
        %    obj = CANDIDATESETSFILE(DATA_DIR, LUT) creates an
        %    object representing a candidate set that is stored as a
        %    csv-file.
        %
        %    NOTATION:
        %       l, scalar: Number of test / independent examples
        %       m, scalar: Number of different candidate sets.
        %       d, scalar: Number of molecular features, e.g. fingerprints,
        %
        %    INPUTS:
        %       data_dir, string
        %               Directory containing the sets of candiate list
        %               csv-files associated with the test or independent
        %               set under investigation.
        %
        %       lut, string array, of length l
        %               Look-Up-Table, mapping each test / independent
        %               examples i in {1,...,l} to one of the m candidate
        %               sets.
        %
        %       selec_feature, logical array, of length d, 'true' indices 
        %               represent the used features. The dimension of the 
        %               candidate feature representation, e.g.,
        %               fingerprints is determined from this array.
        %
        %    OUTPUT: 
        %       obj     A valid object of class CANDIDATESETSFILE.
            if (nargin < 2)
                error('CandidateSetsFiles_2:CandidateSetsFiles:InvalidArgument', 'Not enough input arguments.\n');
            end % if
            
            if (~ exist(data_dir, 'dir')) 
                error('CandidateSetsFiles_2:CandidateSetsFiles:NoSuchFileOrDirectory', ...
                    'Candidate directory "%s" does not exist.', data_dir);
            else
                obj.data_dir_ = data_dir;
            end % if
            
            obj.lut_ = lut;
            obj.selec_feature_ = selec_feature;
            obj.n_features_ = length(obj.selec_feature_);
        end % function
        
        function lhs = getCandidateSet (rhs, id, only_selection, fieldname)
        %% GETCANDIDATESET returns the candidate set for a given example
        %    candSet = GETCANDIDATESET (OBJ, ID) returns the candidate set for 
        %    the example with the given ID.
        %
        %    candSet = GETCANDIDATESET (OBJ, ID, ONLYSELECTION) returns the 
        %    candidate set for the example with the given ID. If 
        %    ONLYSELECTION is 'true', only those elements belonging to 
        %    the currently in the object stored selection are returned. The 
        %    fields storing the representations of the candidates (e.g. 
        %    feature vectors) and the identifier (e.g the inchi) are 
        %    updated according to the selection. The field storing the 
        %    number of candidates in the desired set is updated as well and 
        %    will be equal to the number of selected candidates.
        %
        %    DATA = GETCANDIDATESET (OBJ, ID, ONLYSELECTION, FIELDNAME) returns
        %    the desired field of the candidate set for the example with
        %    the given ID. If ONLYSELECTION is 'true', only the selected 
        %    elements are returned (compare description of the previous call). 
        %
        %    INPUTS:
        %       obj, class-object CandidateSetsFile          
        %
        %       id, integer OR string
        %           Index of the test / independent example for which the
        %           candidate the set information are desired.
        %
        %           OR
        %
        %           ID, e.g. SM854852__01__C19H20ClNO4 of the test /
        %           independent example.
        %       onlySelection, logical (default = false)      
        %           Indicating whether only information from the selected
        %           elements should be returned. 
        %           NOTE: Currently not supported and set to 'false'.
        %       fieldname, string (default = '')
        %           Information about the candidates that should be
        %           returned, e.g. 'fps', 'id', ... 
        %    OUTPUT:
        %       candSet       Structure storing the candidate set
        %                     information.
        %       DATA          Content of the desired field.
        %       
        %       or
        %
        %       NaN           Not a number is returned if no candidate set
        %                     is available for the provided id.
        %       {}            An empty cell is returned if isempty(id) is 
        %                     true.
        %               
        %    EXAMPLE:
        %       candSet_id = Y_c.getCandidateSet(id); 
        %       fp         = Y_c.getCandidateSet(id, 0, 'data'); 
            % Check input
            if (isstring(id) || ischar(id)) 
                id = find(cellfun(@(c) strcmp(c, id), rhs.lut_));
            end % if
            if (isempty(id)); lhs = {}; return; end % if
            if (numel(id) > 1)
                error ('CandidateSetsFile_2:getCandidateSet:InvalidArgument', 'You can only access one element at a time.');
            end % if
            if (id < 1) || (id > numel(rhs.lut_))
                errorStr = 'Index out of bounds';
                if (isempty(rhs.lut_))
                    errorStr = strcat(errorStr, ': LUT is empty.');
                end % if
                error('CandidateSetsFile_2:getCandidateSet:OutOfRange', errorStr);
            end % if 
            
            if (nargin < 3) % NOTE: 'rhs' is a paramter
                only_selection = false;
            end % if
            
            if (only_selection)
                warning('CandidateSetsFile_2:getCandidateSet:UnsupportedFeature', ...
                    'Currently selection of candidates is not supported when loaded from files. "onlySelection" will be set to "false"');
                only_selection = false;
            end % if
            
            if (nargin < 4)
                fieldname = '';
            else
                if (~ ismember(fieldname, fieldnames(CandidateSetsFile_2.candidate_set_fieldnames)))
                    error('CandidateSetFile_2:getCandidateSet:InvalidArgument', ...
                        'Fieldname "%s" does not exist.', fieldname);
                end % if
            end % if
            
            cand_fn = fullfile(rhs.data_dir_, strcat(rhs.lut_(id), ".csv"));
            if (~ exist(cand_fn, 'file'))
                error('CandidateSetsFile_2:getCandidateSet:NoSuchFileOrDirectory', ...
                    'Candidate set file "%s" does not exist.', cand_fn);
            end % if
            
            cand_fid = fopen(cand_fn, 'r');
            if (cand_fid == -1)
                error('CandidateSetsFile_2:getCandidateSet:IOError', 'Cannot open "%s",', cand_fn);
            end % if
            
            if ~ isempty(fieldname)
                % Read individual candidate set information
                switch (fieldname)
                    case 'num'
                        % Number of molecular candidates: We count the
                        % lines in the candidate set file.
                        [~, result] = system(strcat("wc -l <", cand_fn));  % Only works on UNIX
                        num = str2double(result);
                          
%                         num = 0;
%                         while ischar(fgetl(cand_fid))
%                             num = num + 1;
%                         end % while
                        
                        lhs = num;
                    case 'id'
                        % Read only first column, containing the candidate
                        % identifiers.
                        cand_ids = cell(20000, 1);  % Pre-allocate large cell array
                        ii = 1;
                        line = fgetl(cand_fid);
                        while ischar(line)
                            cand_ids{ii} = line(1:14);  % Length of InChIKey1 is 14 characters
                            line = fgetl(cand_fid);
                            ii = ii + 1;
                        end % while
                        cand_ids(cellfun(@isempty, cand_ids)) = [];  % Remove empty cells
                        
                        lhs = cand_ids;
                    case 'data'
                        % Read the information about the features.
                        % These are in all the remaining columns after the
                        % first one. In the case of fingerprints: In
                        % the file might be different number of entries in each row, 
                        % as only active ones fingerprints are included.
                        data = false(rhs.n_features_, 20000);
                        ii = 1;
                        line = fgetl(cand_fid);
                        while ischar(line)
                            data(sscanf(line(16:end), "%d") + 1, ii) = true;
                            line = fgetl(cand_fid);
                            ii = ii + 1;
                        end % while
                        data(:, ii:end) = [];  % Remove empty columns
                        
                        lhs = data;
                end % switch
            else
                error('CandidateSetsFile_2:getCandidateSet:UnsupportedFeature', ...
                    'Currently only specific information from the candidate set can be returned.')
            end % if     
            
            fclose(cand_fid);
        end % function
        
        function numberOfExamples = getNumberOfExamples (obj)
        %% GETNUMBEROFEXAMPLES number of examples associated with the candidate sets    
        %    numberOfExamples = GETNUMBEROFEXAMPLES (OBJ) returns the
        %    number of examples associated with the candidate sets stored
        %    in the OBJ.
        %
        %    OUTPUTS:
        %       numberOfExamples        Intenger 
            numberOfExamples = numel (obj.lut_);
        end % if
        
        function lhs = getSubset (rhs, ids)
        %% GETSUBSET returns an obj. CanidateSets with updated Look-Up-Table and selection
        %    obj = GETSUBSET (RHS, IDS) returns an object CandidateSets,
        %    whereby OBJ's look-up-table contains only the elements
        %    belonging to the IDS. Furthermore, the SELEC vector is updated
        %    and also only contains the elements at position IDS. 
        %
        %    INPUTS:
        %       rhs         Object of class CandidateSets
        %       ids         Indices of the example which should be part of
        %                   the selected subset.
        %
        %    OUTPUT:
        %       lhs         Object of class CandidateSets with an updated
        %                   look-up-table and selection vector.
        %
        %    EXAMPLE:
        %       Y_c = CandidateSets (DataHandle (cand), mf_corres); 
        %       Y_c_train = Y_c.getSubset (idx_train)
        %
        %    TODO: 
        %       * Allow also logical indices to get a subset of the
        %         candidates.
            if (isstring(ids) || ischar(ids)) 
                ids = find(cellfun(@(c) strcmp(c, ids), rhs.lut_));
                if (isempty(ids))
                    error('CandidateSetsFile_2:getSubset:OutOfRange', ...
                        'Cannot find ids in lut.');
                end % if
            end % if
            
            if (any (ids < 1) || any (ids > numel (rhs.lut_)))
                errorStr = 'Index out of bounds';
                if (isempty (rhs.lut_))
                    errorStr = strcat (errorStr, ': LUT is empty.');
                end % if
                error ('CandidateSetsFile_2:getSubset:OutOfRange', errorStr);
            end % if
            
            lhs = CandidateSetsFile_2 (rhs.data_dir_, rhs.lut_(ids), rhs.selec_feature_);
        end % function    
        
        function selec = findExampleInCandidateSet (obj, idx, identifier)
        %% FINDEXAMPLEINCANDIDATESET returns a selection only of the candidate with matching identifier
        %    The returned selection always only contains a single true
        %    value. This means that if several identifier in the candidate
        %    set match the provided IDENTIFIER only the first index is 
        %    returned.
        %
        %    selec = FINDEXAMPLEINCANDIDATESET (OBJ, IDX, IDENTIFIER)
        %    Returns a selection with a single true value. The position of
        %    the true value corresponds to the position of the IDENTIFIER
        %    in the OBJ's candidate set for example IDX. 
        %
        %    INPUTS:
        %       obj             An object of class CANDIDATESETS.
        %       idx             Index of the example it's candidate set
        %                       should be searched for the IDENTIFIER.
        %       identifier      string containing the identifier to search 
        %                       for.
        %
        %    OUTPUTS:
        %       selec           Logical vector of the dimension (1 x n_idx).
        %                       The vector contains only one non-zero
        %                       entry, which corresponds to the _first_
        %                       occurence of the IDENTIFIER in the vector
        %                       of identifiers associated with the
        %                       candidate set of example IDX.
            identifiersOfCandidates = obj.getCandidateSet (idx, false, 'id');
            
            selec = strcmp (identifiersOfCandidates, identifier);
            if (~ any (selec))
                warning('Identifier "%s" not found in candidate set "%s".', identifier, ...
                    obj.lut_{idx});
                return;
            end % if
            
            if (sum (selec) > 1)
                warning('The identifier in the candidate set corresponding to example %d are not all unique.', ...
                    idx);
            end % if
                
            identifierIdx = find(selec, 1); % find (X, K) ... 'K' first
            selec = false(size(selec));
            selec(identifierIdx) = true;

            assert(sum(selec) < 2, 'Upps?! There should be only one "true" in the selection.');
        end % if
    end % public methods
end % class
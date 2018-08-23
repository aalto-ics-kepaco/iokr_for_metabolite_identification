classdef CandidateSetsFile < handle
    properties(SetAccess = private, GetAccess = public)
        data_dir_;
        lut_;
        selec_feature_;
    end % properties: write private, read public
    
    properties(Constant, Access = private)
        candidate_set_fieldnames = struct( ...
            'data', 'fingerprints', ...        
            'id', 'inchi2D', ...
            'num', 'numberOfCandidates');
    end % properties: constant, private
    
    properties(Access = private)
        csv_options = struct('delimiter', '\t', 'headerlines', 1, 'formatstr', '%s%s%s');
    end % properties: private
    
    methods(Access = public)
        function obj = CandidateSetsFile(data_dir, lut, selec_feature, csv_options)
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
        %       lut, string array, (1 x l)
        %               Look-Up-Table, mapping each test / independent
        %               examples i in {1,...,l} to one of the m candidate
        %               sets.
        %
        %       selec_feature, logical array, (1 x d) OR string
        %               The 'true' indices represent the used features.
        %
        %               OR
        %
        %               String 'ALL' if all features should be selected.
        %
        %    OUTPUT: 
        %       obj     A valid object of class CANDIDATESETSFILE.
            if (nargin < 2)
                error('CandidateSetsFiles:CandidateSetsFiles:InvalidArgument', 'Not enough input arguments.\n');
            end % if
            
            if (~ exist(data_dir, 'dir')) 
                error('CandidateSetsFiles:CandidateSetsFiles:NoSuchFileOrDirectory', ...
                    'Candidate directory "%s" does not exist.', data_dir);
            else
                obj.data_dir_ = data_dir;
            end % if
            
            obj.lut_ = lut;
            
            if (nargin > 2)
                obj.selec_feature_ = selec_feature;
            else
                obj.selec_feature_ = 'ALL';
            end % if
            
            if (nargin > 3)
                obj.csv_options = csv_options;
            end % if
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
                error ('CandidateSets:getCandidateSet:InvalidArgument', 'You can only access one element at a time.');
            end % if
            if (id < 1) || (id > numel(rhs.lut_))
                errorStr = 'Index out of bounds';
                if (isempty(rhs.lut_))
                    errorStr = strcat(errorStr, ': LUT is empty.');
                end % if
                error('CandidateSets:getCandidateSet:OutOfRange', errorStr);
            end % if 
            
            if (nargin < 3) % NOTE: 'rhs' is a paramter
                only_selection = false;
            end % if
            
            if (only_selection)
                warning('CandidateSetsFile:getCandidateSet:UnsupportedFeature', ...
                    'Currently selection of candidates is not supported when loaded from files. "onlySelection" will be set to "false"');
                only_selection = false;
            end % if
            
            if (nargin < 4)
                fieldname = '';
            else
                if (~ ismember(fieldname, fieldnames(CandidateSetsFile.candidate_set_fieldnames)))
                    error('CandidateSet:getCandidateSet:InvalidArgument', ...
                        'Fieldname "%s" does not exist.', fieldname);
                end % if
            end % if
            
            cand_fn = rhs.data_dir_ + "/" + rhs.lut_(id) + ".fpt";
            if (~ exist(cand_fn, 'file'))
                error('CandidateSet:getCandidateSet:NoSuchFileOrDirectory', ...
                    'Candidate set file "%s" does not exist.', cand_fn);
            end % if
            
            cand_fid = fopen(cand_fn, 'r');
            if (~ cand_fid)
                error('CandidateSet:getCandidateSet:IOError', 'Cannot open "%s",', cand_fn);
            end % if
            
            % Read header of the candidate file to determine the relevant
            % columns 
            header = textscan(cand_fid, rhs.csv_options.formatstr, 1, 'Delimiter', rhs.csv_options.delimiter);
            
            % Get the column indices for the ids, e.g. inchi, and the data,
            % e.g. fingerprints.
            id_col = cellfun(@(x) strcmp(x, CandidateSetsFile.candidate_set_fieldnames.id), header);
            data_col = cellfun(@(x) strcmp(x, CandidateSetsFile.candidate_set_fieldnames.data), header);
            
            % Read in the whole matrix. 
            % TODO: If only the ids are required, we do not need to read
            % everything.
            cand_data = textscan(cand_fid, rhs.csv_options.formatstr, 'Delimiter', rhs.csv_options.delimiter);
            fclose(cand_fid);
            
            % Transform the fingerprint into sparse matrix
            representation = sparse(cell2mat(cellfun(@(c) c - '0', cand_data{data_col}, 'uniformoutput', false)))';
            
            % Get mask for the fingerprints
            d = size(representation, 1);
            if (strcmpi(rhs.selec_feature_, 'ALL'))
                selec_feature = true(d, 1);
            else
                if (numel(rhs.selec_feature_) ~= d) 
                    error('CandidateSet:getCandidateSet:InvalidArgument', ...
                        'Dimension of the feature selection mask does not math the feature dimension' + ...
                        '(%d vs. %d)', numel(rhs.selec_feature_), d);
                end % if
                selec_feature = rhs.selec_feature_;
            end % if
            
            if (isempty(fieldname)) 
                lhs = struct();
                lhs.data = representation(selec_feature, :);
                lhs.id = cand_data{id_col};
                lhs.num = numel(cand_data{id_col});
            else
                switch (fieldname)
                    case 'id'
                        lhs = cand_data{id_col};
                    case 'data'
                        lhs = representation(selec_feature, :);
                    case 'num'
                        lhs = numel(cand_data{id_col});
                    otherwise
                        error('CandidateSet:getCandidateSet:InvalidArgument', 'No field with name "%s".', ...
                            fieldname)
                end % switch
            end % if 
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
            identifiersOfCandidates = obj.getCandidateSet ( ...
                idx, 0, CandidateSets.candidateSetFieldnames.id);
            
            selec = strcmp (identifiersOfCandidates, identifier);
            if (~ any (selec))
                warning('Identifier not found in candidate set.');        
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
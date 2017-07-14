%% CANDIDATESETS
classdef CandidateSets < handle
    properties (SetAccess = private, GetAccess = public)
        data_handle_;
        lut_;
        selec_;
    end % properties: write private, read public
    
    properties (Constant, Access = private)
        candidateSetFieldnames = struct ( ...
            'representation',   'data', ...        
            'identifier',       'id', ...
            'numberOfElements', 'num');
    end % properties: 
    
    methods (Static, Access = public)      
        function [isValid, errorStr] = validateCandidateSetDataStructure (data_handle)
        %% VALIDATECANDIDATESETDATASTRUCTURE validates structure of candidate data
            isValid = true;
            errorStr = '';
            
            % Check the data handle properties
            if (isempty (data_handle))
                errorStr = 'Candidate set data_handle must not be empty.';
                
                isValid = false; 
                return;
            end % if
            if (~ isa (data_handle, 'DataHandle'))
                errorStr = 'Candidate set data_handle must be a DataHandle.';
                
                isValid = false;
                return;
            end % if
            
            % Check the struct properties
            if (~ isstruct (data_handle.data_))
                errorStr = sprintf ('Candidate set data must be stored in a struct-vector. (provided class = %s)', ...
                    class (data_handle.data_));
                
                isValid = false; 
                return;
            end % if
            if (all (size (data_handle.data_) > 1))
                errorStr = sprintf ('Candidate set data must be stored either as row or column vector. (provided size = %dx%d)', ...
                    size (data_handle.data_, 1), size (data_handle.data_, 2));
                
                isValid = false;
                return; 
            end % if
            if (isempty (data_handle.data_)) % isempty (struct ([])) --> 1
                errorStr = 'Candidate set data must contain at least one candidate set.';
                
                isValid = false; 
                return;
            end % if
            tmp = {CandidateSets.candidateSetFieldnames.representation, ...
                CandidateSets.candidateSetFieldnames.identifier, ...
                CandidateSets.candidateSetFieldnames.numberOfElements};
            if (~ all (ismember (tmp, fieldnames (data_handle.data_))))
                errorStr = 'Candidate set structs must at least contain certain fields (see implementation of CandidateSets)';
                
                isValid = false; 
                return;
            end % if
            if (any (arrayfun (@(a) any ([ ...
                    isempty(a.(CandidateSets.candidateSetFieldnames.representation)), ...
                    isempty(a.(CandidateSets.candidateSetFieldnames.identifier)), ...
                    a.(CandidateSets.candidateSetFieldnames.numberOfElements) == 0]), data_handle.data_)))
                errorStr = 'There should not be any empty candidate set.';
                
                isValid = false;
                return;
            end % if
            if (~ all (arrayfun (...
                    @(a) size (a.(CandidateSets.candidateSetFieldnames.representation), 2) == size (a.(CandidateSets.candidateSetFieldnames.identifier), 2), ...
                    data_handle.data_)))
                errorStr = 'The number of candidate-representations and the number of candidate-identifiers must be the same.';
                
                isValid = false;
                return;
            end % if
            if (~ all (arrayfun (...
                    @(a) size (a.(CandidateSets.candidateSetFieldnames.representation), 2) == a.(CandidateSets.candidateSetFieldnames.numberOfElements), ...
                    data_handle.data_)))
                errorStr = 'The number of candidates must be the same as the number of candidate-representations and candidate-identifiers.';
                
                isValid = false;
                return;
            end % if
            
        end % function 
        
        function [isValid, errorStr] = validateLUT (data_handle, lut)
        %% VALIDATELUT validates look-up-table structure
            isValid = true;
            errorStr = '';
            
            if (size (lut, 2) > 1)
                errorStr = 'Look-Up-Table must be a column-vector.';
                
                isValid = false; 
                return;
            end % if
            if (numel (data_handle.data_) < max (lut))             
                errorStr = 'Look-up-table does contain candidate set indices which does not exists in the data.';
                
                isValid = false;
                return;
            end % if            
        end % function
        
        function [isValid, errorStr] = validateSelection (data_handle, lut, selec)
        %% VALIDATESELECTION validates the selection structure
            isValid = true;
            errorStr = '';
            
            if (~ iscell(selec) || all (size (selec) > 1))
                errorStr = 'Selections must be provided as cell-vector.';
                
                isValid = false;
                return;
            end % if
            if (~ all (size (lut) == size (selec)))
                errorStr = sprintf (...
                    'Dimensions of Look-Up-Table (%d x %d) and selection (%d x %d) must match.', ...
                    size (lut, 1), size (lut, 2), ...
                    size (selec, 1), size (selec, 2));
                
                isValid = false;
                return;
            end % if              
            for i = 1:numel (selec)
                if (size (selec{i}, 1) > 1)
                    errorStr = 'Indices of the selections must be stored in a row-vector.';

                    isValid = false;
                    return;
                end % if
                if (any (isnan (selec{i})) && (numel (selec{i}) > 1))
                    errorStr = 'Elements in selec which correspond to a NaN entry in LUT should be single NaN values.';
                    
                    isValid = false;
                    return;
                end % if
                % TODO: Seems that this does not hold for the metabolite
                % identification data-set.
%                 if (isnan (lut(i)) ~= any (isnan (selec{i})))
%                     errorStr = 'The rows of the "NaN"-entries of lut and selec must match.';
%                     
%                     isValid = false;
%                     return;
%                 end % if
                if (~ any (isnan (selec{i})) & ~ islogical (selec{i}))
                    errorStr = 'The selection vector must be of class logical or NaN.';
                    
                    isValid = false;
                    return; 
                end % if
            end % for
            for i = 1:numel (lut)
                if (~ isnan (lut(i)))
                    if (numel (selec{i}) ~= data_handle.data_(lut(i)).(CandidateSets.candidateSetFieldnames.numberOfElements))
                        errorStr = sprintf ( ...
                            'Length of selec{%d} (=%d) must match the number of candidates in set %d (=%d).', ...
                            i, numel (selec{i}), ...
                            i, data_handle.data_(lut(i)).(CandidateSets.candidateSetFieldnames.numberOfElements));
                        
                        isValid = false;
                        return;
                    end % if
                    if (numel (selec{i}) ~= size (data_handle.data_(lut(i)).(CandidateSets.candidateSetFieldnames.representation), 2))
                        errorStr = sprintf ( ...
                            'Length of selec{%d} (=%d) must match the present representations in candidate set %d (=%d).', ...
                            i, numel (selec{i}), ...
                            i, size (data_handle.data_(lut(i)).(CandidateSets.candidateSetFieldnames.representation), 2));
                        
                        isValid = false;
                        return;
                    end % if
                    if (numel (selec{i}) ~= numel (data_handle.data_(lut(i)).(CandidateSets.candidateSetFieldnames.identifier)))
                        errorStr = sprintf ( ...
                            'Length of selec{%d} (=%d) must match the present indetifiers in candidate set %d (=%d).', ...
                            i, numel (selec{i}), ...
                            i, numel (data_handle.data_(lut(i)).(CandidateSets.candidateSetFieldnames.identifier)));
                        
                        isValid = false;
                        return;
                    end % if
                end %if
            end % for
        end % function
    end % methods
    
    methods (Access = public)
        function obj = CandidateSets (data_handle, lut, selec)
        %% CANDIDATESETS Constructor of the class
        %    obj = CANDIDATESETS (PTRTOCANDIDATEDATA, LUT) creates an object
        %    representing a candidate set. A pointer to the data of the 
        %    candidate sets (e.g. representation, identifier, ...) must be
        %    provided to the function. Furthermore the Look-Up-Table
        %    containing the candidate set id for each example is must be
        %    provided. As no selection for the candidates is provided, all
        %    candidate are selected.
        %    
        %    obj = CANDIDATESETS (PTRTOCANDIATEDATA, LUT, SELECTION) creates 
        %    an object representing a candidate set. The SELECTION defines
        %    for each example the selection of the candidates of the
        %    associated candidate set. 
        %
        %    INPUTS:
        %       ptrToCandidateData 
        %               DATAHANLDE storing the (m x 1) dimensional struct
        %               array containing the candidate sets as structures
        %               of the following shape:
        %               data(i) = struct ('data', DATA, 'id', IDS, ...
        %                   'num', NUMBEROFELEMENTS)
        %       lut     The Look-Up-Table is a (l x 1) dimensional column
        %               vector. Each entry in this vector associates an
        %               example with candidate set. The elements of this
        %               vector are from the range {1, ..., m}. 
        %       selection
        %               (l x 1) dimensional cell array. Each cell i contains
        %               a logical vector of length m(i). This logical
        %               vector indicates whether a certain candidate is in
        %               the selection or not.
        %
        %    OUTPUT: 
        %       obj     A valid object of class CANDIDATESETS. 
        %
        %    See also DATAHANDLE
            if (nargin < 1)
                error ('CandidateSets:CandidateSets:InvalidArgument', ...
                    'Not enough input arguments.\n');
            end % if
            if (nargin < 2)
                lut = [];
            end % if
        
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure (data_handle);
            if (~ isValid)
                error ('CandidateSets:CandidateSets:InvalidArgument', errorStr);
            end % if
            
            obj.data_handle_ = data_handle;
            
            [isValid, errorStr] = CandidateSets.validateLUT (obj.data_handle_, lut);
            if (~ isValid)
                error ('CandidateSets:CandidateSets:InvalidArgument', errorStr);
            end % if
            
            obj.lut_ = lut;
            
            % Selection is provided            
            if (nargin > 2)
                [isValid, errorStr] = CandidateSets.validateSelection (obj.data_handle_, obj.lut_, selec);
                if (~ isValid)
                    error ('CandidateSets:CandidateSets:InvalidArgument', errorStr);
                end % if
                            
                obj.selec_ = selec;
            else
                obj.resetSelectionToAllCandidates();
            end % if
        end % function: constructor
        
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
            if (islogical (ids))
                if (numel (ids) ~= numel (rhs.lut_) ...
                    || numel (ids) ~= numel (rhs.selec_))
                    error ('CandidateSets:getSubset:OutOfRange', ...
                        'Logical index-vector does not match the dimension of "lut" and "selec".');
                end % if
            else
                if (any (ids < 1) || any (ids > numel (rhs.lut_)))
                    errorStr = 'Index out of bounds';
                    if (isempty (rhs.lut_))
                        errorStr = strcat (errorStr, ': LUT is empty.');
                    end % if
                    error ('CandidateSets:getSubset:OutOfRange', errorStr);
                end % if
            end % if
            
            lhs = CandidateSets (rhs.data_handle_, rhs.lut_(ids), rhs.selec_(ids));
        end % function    
        
        function lhs = getCandidateSet (rhs, id, onlySelection, fieldname)
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
        %       obj           Is a valid CandidateSets object.
        %       id            Integer id of the example for which 
        %                     information is desired.
        %       onlySelection Binary indicating whether only information
        %                     from the selected elements should be
        %                     returned.                   [default = false]
        %       fieldname     String naming the fieldname which should be
        %                     returned.                      [default = '']
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
        %       candSet_id = Y_c.getElement(id); 
        %       fp         = Y_c.getElement(id, 0, 'data'); 
        %
        %    See also GETCANDIDATESELECTION     
        
            % Check input
            if (isempty (id));         lhs = {};  return; end % if
            if (isnan (rhs.lut_(id))); lhs = NaN; return; end % if
            
            if (numel (id) > 1)
                error ('CandidateSets:getCandidateSet:InvalidArgument', ...
                    'You can only access one element at a time.');
            end % if
            if (id < 1) || (id > numel (rhs.lut_))
                errorStr = 'Index out of bounds';
                if (isempty (rhs.lut_))
                    errorStr = strcat (errorStr, ': LUT is empty.');
                end % if
                error ('CandidateSets:getCandidateSet:OutOfRange', errorStr);
            end % if
            
            if (nargin < 3) % NOTE: 'rhs' is a paramter
                onlySelection = false;
            end % if
            if (nargin < 4)
                fieldname = '';
            else
                if (~ ismember (fieldname, fieldnames (rhs.data_handle_.data_)))
                    error ('CandidateSet:getCandidateSet:InvalidArgument', ...
                        'Fieldname "%s" does not exist.', fieldname);
                end % if
            end % if
            
            % Access element            
            switch (onlySelection) 
                case false
                    if (isempty (fieldname))
                        lhs = rhs.data_handle_.data_(rhs.lut_(id));
                    else
                        lhs = rhs.data_handle_.data_(rhs.lut_(id)).(fieldname);
                    end % if
                case true
                    if (isempty (fieldname))
                        % Copy the tree main fields into the output.
                        lhs = struct ();
                        lhs.(CandidateSets.candidateSetFieldnames.representation) = ...
                            rhs.data_handle_.data_(rhs.lut_(id)).(CandidateSets.candidateSetFieldnames.representation)(:, rhs.selec_{id});
                        lhs.(CandidateSets.candidateSetFieldnames.identifier) = ...
                            rhs.data_handle_.data_(rhs.lut_(id)).(CandidateSets.candidateSetFieldnames.identifier)(rhs.selec_{id});
                        lhs.(CandidateSets.candidateSetFieldnames.numberOfElements) = sum (rhs.selec_{id});
                        
                        % Copy the remaining fields into the output without
                        % any selection.
                        tmp = {CandidateSets.candidateSetFieldnames.representation, ...
                            CandidateSets.candidateSetFieldnames.identifier, ...
                            CandidateSets.candidateSetFieldnames.numberOfElements};
                        remainingFn = setdiff (fieldnames (rhs.data_handle_.data_(rhs.lut_(id))), tmp);
                        for fn = remainingFn'
                            lhs.(fn{1}) = rhs.data_handle_.data_(rhs.lut_(id)).(fn{1});
                        end % for
                    else
                        switch (fieldname)
                            case CandidateSets.candidateSetFieldnames.representation
                                lhs = rhs.data_handle_.data_(rhs.lut_(id)).(fieldname)(:, rhs.selec_{id});
                            case CandidateSets.candidateSetFieldnames.identifier
                                lhs = rhs.data_handle_.data_(rhs.lut_(id)).(fieldname)(rhs.selec_{id});
                            case CandidateSets.candidateSetFieldnames.numberOfElements
                                lhs = sum (rhs.selec_{id});
                            otherwise
                                lhs = rhs.data_handle_.data_(rhs.lut_(id)).(fieldname);
                        end % switch
                    end % if
                otherwise 
                    error ('CandidateSets:getCandidateSet:InvalidArgument', ...
                        '"onlySelection" must be binary.');
            end % switch
            
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
        
        function resetSelectionToAllCandidates (obj)
        %% RESETSELECTIONTOALLCANDIDATES resets the selection of candidates to all
            selec = obj.createSelectionOfAllCandidates();

            [isValid, errorStr] = CandidateSets.validateSelection (obj.data_handle_, obj.lut_, selec);
            assert (isValid, 'CandidateSets:createSelectionOfAllCandidates Ups!?: %s', errorStr);

            obj.selec_ = selec;
        end % function
        
        function setSelectionsOfCandidateSets (obj, selec)
        %% SETSELECTIONSOFCANDIDATESETS sets the SELEC parameter for all candidate sets
            [isValid, errorStr] = CandidateSets.validateSelection ( ...
                obj.data_handle_, obj.lut_, selec);
            if (~ isValid)
                error ('CandidateSets:setSelectionsOfCandidateSets:InvalidArgument', ...
                    'Message: %s', errorStr);
            end % if
            
            obj.selec_ = selec;
        end % function
        
        function selec = getSelectionsOfCandidateSets (obj)
        %% GETSELECTIONSOFCANDIDATESETS returns the SELEC parameter for all candidate sets
            selec = obj.selec_;
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
            identifiersOfCandidates = obj.getCandidateSet ( ...
                idx, 0, CandidateSets.candidateSetFieldnames.identifier);
            
            selec = strcmp (identifiersOfCandidates, identifier);
            if (~ any (selec))
%                 warning ('Identifier not found in candidate set.');
                
                return;
            end % if
            
            if (sum (selec) > 1)
                warning ('The identifier in the candidate set corresponding to example %d are not all unique.', ...
                    idx);
            end % if
                
            identifierIdx = find (selec, 1); % find (X, K) ... 'K' first
            selec = false (size (selec));
            selec(identifierIdx) = true;

            assert (sum (selec) < 2, 'Upps?! There should be only one "true" in the selection.');
        end % if
        
        function selec = createSelectionOfAllCandidates (obj)
        %% CREATESELECTIONOFALLCANDIDATES creates a selection of all candidates
            selec = cell (size (obj.lut_));
                
            for i = 1:numel (obj.lut_)
                if (~ isnan (obj.lut_(i)))
                    selec{i} = true (...
                        1, obj.data_handle_.data_(obj.lut_(i)).(CandidateSets.candidateSetFieldnames.numberOfElements));
                else
                    selec{i} = NaN;
                end % if
            end % for
            
            assert (all (cellfun (@(c) xor (all (isnan (c)), all (c == true)), selec, 'UniformOutput', true)), 'Upps?!');
        end % function
        
        function isEqual = eq (obj1, obj2)
        %% ISEQUAL returns true if two CanidateSets objects are equal.
            isEqual = 1;
            
            % some basics
            isEqual = isEqual & (isempty (obj1) == isempty (obj2));
            if (~ isEqual) ; return ; end % if
            isEqual = isEqual & (strcmp (class (obj1), class (obj2)));
            if (~ isEqual) ; return ; end % if
            
            % dimensions
            isEqual = isEqual & (numel (obj1.lut_) == numel (obj2.lut_));
            if (~ isEqual) ; return ; end % if
            isEqual = isEqual & (numel (obj1.selec_) == numel (obj2.selec_));
            if (~ isEqual) ; return ; end % if
            isEqual = isEqual & (all (cellfun (@(c1, c2) all (size (c1) == size (c2)), obj2.selec_, obj1.selec_)));
            if (~ isEqual) ; return ; end % if
            isEqual = isEqual & (numel (obj1.data_handle_.data_) == numel (obj1.data_handle_.data_));
            if (~ isEqual) ; return ; end % if
            
            % content: lut and selec
            isEqual = isEqual & (all (obj1.lut_ == obj2.lut_));
            if (~ isEqual) ; return ; end % if
            isEqual = isEqual & (all (cellfun (@(c1, c2) all(c1 == c2), obj1.selec_, obj2.selec_)));
            if (~ isEqual) ; return ; end % if
            
            % content: data_handle
            %   fieldnames
            isEqual = isEqual & ...
                (numel (fieldnames (obj1.data_handle_.data_)) == numel (fieldnames (obj2.data_handle_.data_)));
            if (~ isEqual) ; return ; end % if
            isEqual = isEqual & (all (strcmp ( ...
                sort (fieldnames (obj1.data_handle_.data_)), ...
                sort (fieldnames (obj2.data_handle_.data_)))));
            if (~ isEqual) ; return ; end % if
            
            for i = 1:numel (obj1.data_handle_.data_)
                % Compare only the required fields, whether they have
                for fn = fieldnames (CandidateSets.candidateSetFieldnames)'
                    fc = CandidateSets.candidateSetFieldnames.(fn{1});
                    % ... the same class
                    isEqual = isEqual & strcmp ( ...
                        class (obj1.data_handle_.data_(i).(fc)), class (obj2.data_handle_.data_(i).(fc)));
                    if (~ isEqual) ; return ; end % if
                    
                    % ... are storing the same content
                    if (isnumeric (obj1.data_handle_.data_(i).(fc)) || islogical (obj1.data_handle_.data_(i).(fc)))
                        isEqual = isEqual & (all (all (obj1.data_handle_.data_(i).(fc) == obj2.data_handle_.data_(i).(fc))));
                    elseif (iscell (obj1.data_handle_.data_(i).(fc)))
                        if (ischar (obj1.data_handle_.data_(i).(fc){1}))
                            isEqual = isEqual & (all (cellfun (@(c1, c2) strcmp (c1, c2), ...
                                obj1.data_handle_.data_(i).(fc), ...
                                obj2.data_handle_.data_(i).(fc))));
                        else
                            isEqual = isEqual & (all (cellfun (@(c1, c2) c1 == c2, ...
                                obj1.data_handle_.data_(i).(fc), ...
                                obj2.data_handle_.data_(i).(fc))));
                        end % if
                    end % if
                    if (~ isEqual) ; return ; end % if
                end % for
            end % for
        end % function
    end % methods: public
end % class

%% How to use this stuff:
% cand      = load (filename1); 
% mf_corres = load (filename2); 
% 
% Y_c = CandidateSets (DataHandle (cand), mf_corres); 
% Y_c_train = Y_c (1:100)
% 
% clear cand mf_corres; 

%% How to overload an operator
% See also: https://se.mathworks.com/matlabcentral/answers/153131-problem-with-subreferencing-subsref-operator-overload
% function lhs = subsref (rhs, S)
% %             if (any (ids.subs{1} < 1) || any (ids.subs{1} > numel (rhs.lut_)))
% %                 error ('Index out of bounds.');
% %             end % if
%             
%     switch (S.type)
%         case '()'
%             % (ids) returns a new object with a copy of the data
%             % handle, but a reduced look-up-table. 
%             %
%             % EXAMPLE:
%             %   Y_c = CandidateSets (DataHandle (cand), mf_corres); 
%             %   Y_c_train = Y_c (1:100)
%             lhs = CandidateSets (rhs.data_handle_, rhs.lut_(S.subs));
%         case '{}'
%             % {id} returns the candidate set for training example
%             % id. 
%             lhs = rhs.data_handle_.data_{rhs.lut_(S.subs)};
%         case '.'
%             lhs = rhs.(S.subs);
%     end % switch            
% end % function: operator(), operator{} and operator.
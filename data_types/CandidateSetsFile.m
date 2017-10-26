classdef CandidateSetsFile < handle
    methods (Access = public)
        function obj = CandidateSetsFile (data_dir, idx2setdesc, selec, selec_feature)
        %% CANDIDATESETSFILE Class-constructor
        %    obj = CANDIDATESETSFILE (DATA_DIR, IDX2SETDESC) creates an
        %    object representing a candidate set. 
        %
        %    INPUTS:
        %       data_dir
        %               Directory containg a set of csv-files. Those sets which
        %               correspond to the set of candidates
        
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
        %    obj = CANDIDATESETS (PTRTOCANDIATEDATA, LUT, SELECTION, DATA_SELECTION) 
        %    creates an object representing a candidate set. The DATA_SELECTION 
        %    defines the subset of the feature representation returned on 
        %    access: data(DATA_SELECTION, :) returns only a subset of the
        %    candidate representation.
        %
        %    INPUTS:
        %       ptrToCandidateData 
        %               DATAHANLDE storing the (m x 1) dimensional struct
        %               array containing the candidate sets as structures
        %               of the following shape:
        %               data(i) = struct ('data', DATA, 'id', IDS, ...
        %                   'num', NUMBEROFELEMENTS)
        %       lut     The Look-Up-Table is a l - dimensional vector. 
        %               Each entry in this vector associates an
        %               example with a candidate set. The elements of this
        %               vector are from the range {1, ..., m}. 
        %               lut : {1, ..., l} --> {1, ..., m}
        %       selection
        %               (l x 1) dimensional cell array. Each cell i contains
        %               a logical vector of length m(i). This logical
        %               vector indicates whether a certain candidate is in
        %               the selection or not.
        %   
        %               OR
        %
        %               String 'ALL' if all candidate should be selected.
        %               This option is useful if no candidate but a feature
        %               selection is desired.
        %
        %       feature_selection
        %               d - dimensional logical vector. The 'true'
        %               indices represent the used features.
        %
        %               OR
        %
        %               String 'ALL' if all features should be selected.
        %
        %    OUTPUT: 
        %       obj     A valid object of class CANDIDATESETS. 
        %
        %    See also DATAHANDLE
        end % function
    end % public methods
end % class
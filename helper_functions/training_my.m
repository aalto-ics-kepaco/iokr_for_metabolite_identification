%% TRAINING_MY behaves like the cvpartitions/training function
%    The purpose of this function is to provide a joint interface for the
%    cv structures I introduced.
%
%    FOLDSELEC = TRAINING (CV, I) returns a binary vector with the length
%    equal to the number of examples in fold I defined in CV.
%
%    INPUTS:
%       cv          Either a cvpartition object or a struct of the following
%                   structure: 
%                            * training:    (nTrain x nOuterFolds) logical
%                                           matrix. (:, i) corresponds to
%                                           the selection of the training
%                                           examples for outer fold i.
%                            * test:        (nTrain x nOuterFolds) logical
%                                           matrix. (:, i) corresponds to
%                                           the selection of the test
%                                           examples for outer fold i.
%
%    OUTPUTS:
%       foldSelec   (nTrain_i x 1)-dimensional binary vector. The logical
%                   value of each element j indicates whether example j is
%                   part of the training (true) or not (false).
%
%    See also GETCVINDICES, CVPARTITION and TEST_MY.
function foldSelec = training_my (cv, i)
    if (nargin < 2)
        error ('training_my:InvalidArgument', 'Not enough input arguments.');
    end %if

    switch (class (cv))
        case 'cvpartition'
            foldSelec = cv.training (i);
        case 'struct'
            if (~ ismember ('training', fieldnames (cv)))
                error ('training:InvalidArgument', ...
                    'CV structure must contain a field "training".');
            end % if
            
            foldSelec = cv.training (:, i);
        otherwise
            error ('training:InvalidArgument', ...
                'Class of CV must be CVPARTITION or STRUCT.');
    end % switch
end % function
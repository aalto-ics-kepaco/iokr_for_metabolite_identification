function foldSelec = test_my (cv, i)
%% TEST_MY behaves like the cvpartitions/test function
%    The purpose of this function is to provide a joint interface for the
%    cv structures I introduced. 
%
%    FOLDSELEC = TEST_MY (CV, I) returns a binary vector with the length
%    equal to the number of examples in fold I defined in CV.
%
%    INPUTS:
%       cv          Either a cvpartition object or a struct of the following
%                   structure: 
%                            * test:        (nTrain x nOuterFolds) logical
%                                           matrix. (:, i) corresponds to
%                                           the selection of the test
%                                           examples for outer fold i.
%
%    OUTPUTS:
%       foldSelec   (nTrain_i x 1)-dimensional binary vector. The logical
%                   value of each element j indicates whether example j is
%                   part of the test (true) or not (false).
%
%    See also GETCVINDICES, CVPARTITION and TRAINING_MY.
    if (nargin < 2)
        error ('test_my:InvalidArgument', 'Not enough input arguments.');
    end %if

    switch (class (cv))
        case 'cvpartition'
            foldSelec = cv.test (i);
        case 'struct'
            if (~ ismember ('test', fieldnames (cv)))
                error ('training:InvalidArgument', ...
                    'CV structure must contain a field "test".');
            end % if
            
            foldSelec = cv.test (:, i);
        otherwise
            error ('test_my:InvalidArgument', ...
                'Class of CV must be CVPARTITION or STRUCT.');
    end % switch
end % function
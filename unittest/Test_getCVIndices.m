%% TEST_GETCVINDICES methods to tests the "getCVIndices.m" implementation
%
%   See also GETCVINDICES
classdef Test_getCVIndices < matlab.unittest.TestCase   
    methods (Test)
        function testInputOfFixedCVIndices (tc)
        %% TESTINPUTOFFIXEDCVINDICES
        
            % Load some data
            celineWorkDir = ...
                '/m/cs/scratch/kepaco/celine/Metabolites_identification/GNPS/data/';
            ind_fold = load (strcat (celineWorkDir, 'cv_ind.txt'));
            
            cvParam = struct ( ...
                'outer', struct ('type', 'fixed', 'cvInd', ind_fold), ...
                'inner', struct ('nFolds', 10));
            
            cv = getCVIndices (cvParam);
            for i = 1:10
                test_set = test_my (cv.outer, i);
                train_set = training_my (cv.outer, i);
                tc.assertEqual (test_set, ind_fold == i, 'Test selections should be equal.');
                tc.assertEqual (train_set, ind_fold ~= i, 'Train selections should be equal.');
                tc.assertTrue (all (train_set ~= test_set), 'Test and train selection must not overlap');
                
                for j = 1:10
                    test_set_cv = test_my (cv.inner{i}, j);
                    train_set_cv = training_my (cv.inner{i}, j);

                    tc.assertTrue (sum (train_set) == (sum (test_set_cv) + sum (train_set_cv)), ...
                        '# of outer training examples must be equal to the # of inner training and test examples.');
                    tc.assertTrue (all (train_set_cv ~= test_set_cv), ...
                        'Training and test of the inner fold must not overlap.');
                    tc.assertTrue (all (train_set_cv | test_set_cv), ...
                        'Training and test of the inner fold must cover all the outer training examples');
                    tc.assertTrue (all (train_set(train_set) == (train_set_cv | test_set_cv)), ...
                        'All training examples from the outer fold must be used by the inner fold.');
                end % for
            end % for
        end % function
        
        function testOnlyOuterFold (tc)
        %% TESTONLYOUTERFOLD 
        %    If the parameter struct only contains an 'outer' field, the
        %    final cross-validation does not contain information for an
        %    inner cross-validation.
            cvParam = struct ('nObservations', 1000,  'outer', struct ( ...
                'type', 'random', 'nFolds', 10));
            cv = getCVIndices (cvParam);
            
            tc.assertFalse (ismember ('inner', fieldnames (cv)))
            
            cvParam.inner = struct ('nFolds', 10);
            cv = getCVIndices (cvParam);
            
            tc.assertTrue (ismember ('inner', fieldnames (cv)))
            
            clear cvParam cv;
        end % function    
    end % method
end % class
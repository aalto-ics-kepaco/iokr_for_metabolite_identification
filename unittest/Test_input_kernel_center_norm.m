%% TEST_GETCVINDICES methods to tests the "getCVIndices.m" implementation
%
%   See also GETCVINDICES
classdef Test_input_kernel_center_norm < matlab.unittest.TestCase   
    methods (Test)
        function testCenteringAndNormalization (tc)
            tol = 1e-14; 
            
            X = rand (100, 1000) * 100 - 50;
            n = size (X, 2);
            cv = cvpartition (n, 'KFold', 5);
            
            for ii = 1:cv.NumTestSets
                X_train = X(:, cv.training (ii));
                X_test  = X(:, cv.test (ii));
                
                mean_X_train = mean (X_train, 2);
                Phi_X_train  = norma (X_train, mean_X_train, true);
                Phi_X_test   = norma (X_test, mean_X_train, true);
                
                % Reference
                KX_train_ref      = Phi_X_train' * Phi_X_train;
                KX_train_test_ref = Phi_X_train' * Phi_X_test;
                KX_test_ref       = Phi_X_test' * Phi_X_test; 
                
                
                tc.assertTrue (all ((diag (KX_train_ref) - 1) < tol));
                tc.assertTrue (all ((diag (KX_test_ref)  - 1) < tol));
                
                % 
                KX = X' * X;
                [KX_train, KX_train_test, KX_test] = input_kernel_center_norm ( ...
                    KX, cv.training (ii), cv.test (ii), true);
                
                tc.assertTrue (all ((diag (KX_train) - 1) < tol));
                tc.assertTrue (all ((diag (KX_test)  - 1) < tol));
                
                tc.assertTrue (all ((KX_train_ref(:) - KX_train(:)) < tol));
                tc.assertTrue (all ((KX_train_test_ref(:) - KX_train_test(:)) < tol));
                tc.assertTrue (all ((KX_test_ref(:) - KX_test(:)) < tol));
            end % for
        end % function
    end % methods
end % class
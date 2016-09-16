classdef Test_CandidateSets < matlab.unittest.TestCase   
    methods (Test)
        function testValidateCandidateSetDataStructure (tc)
        %% TESTVALIDATECANDIDATESETDATASTRUCTURE test the validation of the candidate set data
        
            % An empty candidate data is not valid
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure([]);
            disp (errorStr);
            tc.assertFalse (isValid);
            
            % Candidate data must be provided as data_handle
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure('not empty');
            disp (errorStr);
            tc.assertFalse (isValid);
            
            % Candidate data must be stored in a struct-array
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle('no struct-array'));
            disp (errorStr);
            tc.assertFalse (isValid);
            
            % Candidate data must contain at least one candidate set
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (struct([])));
            disp (errorStr);
            tc.assertFalse (isValid);
            
            % Candidate data must be either a row or a column vector
            s = struct ([]);
            s(2, 2).a = 1;
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (s));
            disp (errorStr);
            tc.assertFalse (isValid);
                   
            % Candidate sets contained in the data must have the correct
            % structure
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (s));
            disp (errorStr);
            tc.assertFalse (isValid);            
            candValid    = struct ('data', 1, 'id', {'a'}, 'num', 1);
            tc.assertTrue (CandidateSets.validateCandidateSetDataStructure(DataHandle (candValid)));
            tc.assertTrue (CandidateSets.validateCandidateSetDataStructure(DataHandle (vertcat (candValid, candValid))));
            
            % Candidate sets must not contain empty data. We should use 
            % lut(id) = NaN for this.
            candEmpty = struct ('data', zeros(0), 'id', {cell(0)}, 'num', 0);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candEmpty));
            disp (errorStr);
            tc.assertFalse (isValid);
            candEmpty = struct ('data', zeros(1), 'id', {cell(0)}, 'num', 0);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candEmpty));
            disp (errorStr);
            tc.assertFalse (isValid);
            candEmpty = struct ('data', zeros(1), 'id', {cell(1)}, 'num', 0);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candEmpty));
            disp (errorStr);
            tc.assertFalse (isValid);
            candEmpty = struct ('data', zeros(0), 'id', {cell(1)}, 'num', 0);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candEmpty));
            disp (errorStr);
            tc.assertFalse (isValid);
            candEmpty = struct ('data', zeros(0), 'id', {cell(1)}, 'num', 1);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candEmpty));
            disp (errorStr);
            tc.assertFalse (isValid);
            candEmpty = struct ('data', zeros(0), 'id', {cell(0)}, 'num', 1);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candEmpty));
            disp (errorStr);
            tc.assertFalse (isValid);
            candEmpty = struct ('data', zeros(1), 'id', {cell(0)}, 'num', 0);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candEmpty));
            disp (errorStr);
            tc.assertFalse (isValid);
            
            candEmpty = struct ('data', zeros(1), 'id', {cell(1)}, 'num', 1);
            tc.assertTrue (CandidateSets.validateCandidateSetDataStructure(DataHandle (candEmpty)));
            
            % at least 'data', 'id' and 'num' must be there
            candValid = struct ('data', 1, 'id', {'a'}, 'num', 1, 'bla', 99); 
            tc.assertTrue (CandidateSets.validateCandidateSetDataStructure(DataHandle (candValid)));
            
            % lets try a not valid structure
            candNotValid = struct ('data', 0, 'id', {'a'}, 'iq', 120);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candNotValid));
            disp (errorStr);
            tc.assertFalse (isValid);           
            
            candValid = struct ('data', 1, 'id', {'a'}, 'num', 1);
            candSparse = struct ('data', sparse(1), 'id', {'a'}, 'num', 1); % also sparse as data is fine
            tc.assertTrue (CandidateSets.validateCandidateSetDataStructure(DataHandle (vertcat(candSparse, candValid))));      
            tc.assertTrue (CandidateSets.validateCandidateSetDataStructure(DataHandle (vertcat(candSparse))));         
            
            % dimensions / values of 'data', 'id' and 'num' must match
            candValid = struct ('data', ones(1, 2), 'id', {{'a', 'b'}}, 'num', 2);
            tc.assertTrue (CandidateSets.validateCandidateSetDataStructure(DataHandle (candValid)));
            
            candNotValid = struct ('data', ones(1, 3), 'id', {{'a', 'b'}}, 'num', 3);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candNotValid));
            disp (errorStr);
            tc.assertFalse (isValid); 
            
            candNotValid = struct ('data', ones(1, 3), 'id', {{'a', 'b', 'c'}}, 'num', 2);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candNotValid));
            disp (errorStr);
            tc.assertFalse (isValid); 
            
            candNotValid = struct ('data', ones(1, 2), 'id', {{'a', 'b', 'c'}}, 'num', 3);
            [isValid, errorStr] = CandidateSets.validateCandidateSetDataStructure(DataHandle (candNotValid));
            disp (errorStr);
            tc.assertFalse (isValid); 
            
        end % function 
        
        function testValidateSelection (tc)
        %% TESTVALIDATESELECTION test the validation of the selection
            
            % The SELECTION must be provided as cell-vector and its
            % dimensions must match the dimension of the LUT
            candValid = struct ('data', 1, 'id', {'a'}, 'num', 1);
            tc.assertTrue (CandidateSets.validateSelection (DataHandle (candValid), [], {}));
            
            [isValid, errorStr] = CandidateSets.validateSelection (DataHandle (candValid), [], {true(0); true(0)});
            disp (errorStr);
            tc.assertFalse (isValid);
            
            tc.assertTrue (CandidateSets.validateSelection (DataHandle (candValid), [1, 1]', {true(1); true(1)}));
            
            [isValid, errorStr] = CandidateSets.validateSelection (DataHandle (candValid), ones (2, 2), cell(2, 2));
            disp (errorStr);
            tc.assertFalse (isValid);
            
            % The SELECTION must have a NaN iff LUT has NaN for the
            % corresponding example
            candValid = struct ('data', 1, 'id', {'a'}, 'num', 1);
            [isValid, errorStr] = CandidateSets.validateSelection (DataHandle (candValid), [1, NaN]', {true(1); true(0)});
            disp (errorStr);
            tc.assertFalse (isValid);
            
            [isValid, errorStr] = CandidateSets.validateSelection (DataHandle (candValid), [1, NaN]', {true(1), [NaN, NaN]}');
            disp (errorStr);
            tc.assertFalse (isValid);
            
            [isValid, errorStr] = CandidateSets.validateSelection (DataHandle (candValid), [1, NaN]', {true(1), [NaN, true]}');
            disp (errorStr);
            tc.assertFalse (isValid);
            
            tc.assertTrue (CandidateSets.validateSelection (DataHandle (candValid), [NaN, 1]', {NaN, true(1)}'));
            
            % The SELECTION must have the correct length and should be
            % vectors
            s = struct('data', [], 'id', {}, 'num', []); 
            s(1) = struct ('data', 1, 'id', {'a'}, 'num', 1);
            s(2) = struct ('data', sparse(1, 10), 'id', {repmat({'a'}, [1, 10])}, 'num', 10);
            s(3) = struct ('data', ones(1, 7), 'id', {repmat({'a'}, [1, 7])}, 'num', 7);
            tc.assertTrue (CandidateSets.validateSelection (DataHandle (s(1:2)), [1, 2]', {true(1), true(1, 10)}'));
            tc.assertTrue (CandidateSets.validateSelection ( ...
                DataHandle (s), [1, 2, 1, 3]', {true(1), true(1, 10), true(1), true(1, 7)}'));
            
            [isValid, errorStr] = CandidateSets.validateSelection (DataHandle (s(1:2)), [1, 2]', {true(1), true(1, 9)}');
            disp (errorStr);
            tc.assertFalse (isValid);
            [isValid, errorStr] = CandidateSets.validateSelection (DataHandle (s(1:2)), [1, 2]', {1, true(1, 10)}');
            disp (errorStr);
            tc.assertFalse (isValid);
            
            [isValid, errorStr] = CandidateSets.validateSelection (DataHandle (s(1:2)), [1, 2]', {1, true(2, 10)}');
            disp (errorStr);
            tc.assertFalse (isValid);
            
            tc.assertTrue (CandidateSets.validateSelection (DataHandle (s(1:2)), [1, 2, NaN]', {true(1), true(1, 10), NaN}'));
        end % function 
            

        function testConstructor (tc)
        %% TESTCONSTRUCTOR runs test on the correctness of the constructor 
        %    The test performed for the constructor must ensure that the
        %    created object is always a valid object. If no valid object
        %    can be created, the constrcutor should throw an error. 
        
            % Not enough arguments should throw an error
            tc.assertError (@()CandidateSets,             'CandidateSets:CandidateSets:InvalidArgument');
            tc.assertError (@()CandidateSets(DataHandle), 'CandidateSets:CandidateSets:InvalidArgument');
            
            % Empty / not valid arguments for DATA should throw an error
            tc.assertError (@()CandidateSets(DataHandle ([])), 'CandidateSets:CandidateSets:InvalidArgument');
            
            candValid    = struct ('data', 1, 'id', {'a'}, 'num', 1);
            candNotValid = struct ('id', {cell(0)}, 'data', zeros(0), 'num', 0);
            tc.assertInstanceOf (CandidateSets(DataHandle (candValid)), 'CandidateSets');
            tc.assertError (@()CandidateSets(DataHandle (candNotValid)), 'CandidateSets:CandidateSets:InvalidArgument');
            tc.assertError (@()CandidateSets(DataHandle (vertcat(candValid, candNotValid))), 'CandidateSets:CandidateSets:InvalidArgument');
            
            tc.assertInstanceOf (CandidateSets(DataHandle (candValid), []), 'CandidateSets');
            tc.assertInstanceOf (CandidateSets(DataHandle (candValid), [], {}), 'CandidateSets');
     
            candNotValid = struct ('data', 1, 'id', {'a'}, 'bla', 0);
            tc.assertError (@()CandidateSets(DataHandle (candNotValid)), 'CandidateSets:CandidateSets:InvalidArgument');
            
            % Not valid arguments for LUT should throw an error
            candValid = struct ('data', 1, 'id', {'a'}, 'num', 1);
            % LUT is not row vector
            tc.assertError (@()CandidateSets(DataHandle (candValid), [1, 1, 1] ), 'CandidateSets:CandidateSets:InvalidArgument');
            % LUT does containg candidate set indices which do not exist
            tc.assertError (@()CandidateSets(DataHandle (candValid), [1, 2, 1]'), 'CandidateSets:CandidateSets:InvalidArgument');
            
            tc.assertInstanceOf (CandidateSets(DataHandle (candValid), [1, 1, 1]'), 'CandidateSets');
            tc.assertInstanceOf (CandidateSets(DataHandle (vertcat(candValid, candValid)), [1, 2, 1]'), 'CandidateSets');
            
            % Not valid arguments for SELEC should throw an eror
            candValid1 = struct ('data', ones(1, 10), 'id', {repmat({'a'}, [1, 10])}, 'num', 10);
            candValid2 = struct ('id', {repmat({'a'}, [1, 5])}, 'data', sparse(1, 5), 'num', 5);
            candValid3 = struct ('data', 1, 'id', {'a'}, 'num', 1);
            lutValid   = [1, 1, 2, 3]';
            selecValid = {true(1, 10), true(1, 10), true(1, 5), true(1, 1)}';
            
            tc.assertInstanceOf (CandidateSets (DataHandle (vertcat (candValid1, candValid2, candValid3)), lutValid), 'CandidateSets');
            
            candSet1 = CandidateSets ( ...
                DataHandle (vertcat (candValid1, candValid2, candValid3)), ...
                lutValid, selecValid);
            candSet2 = CandidateSets ( ...
                DataHandle (vertcat (candValid1, candValid2, candValid3)), ...
                lutValid);
            tc.assertInstanceOf (candSet1, 'CandidateSets');
            tc.assertInstanceOf (candSet2, 'CandidateSets');
            
            % This might not be a good test: I think the tests should be
            % independent of each other. This means I should be able to
            % test the constructor without assuming, that "eq" is properly
            % implemented. This is because the tests for "eq" require a
            % correct implementation of the constructor. 
%             tc.assertTrue (eq (candSet1, candSet2))   
            
            selecNotValid = {true(1, 10), true(1, 10), true(1, 5), NaN(1, 2)}';
            tc.assertError(@()CandidateSets (DataHandle (vertcat (candValid1, candValid2, candValid3)), lutValid, selecNotValid), ...
                'CandidateSets:CandidateSets:InvalidArgument');
            
            selecNotValid = {true(1, 10), true(1, 10), true(1, 5), true(2, 1)}';
            tc.assertError(@()CandidateSets (DataHandle (vertcat (candValid1, candValid2, candValid3)), lutValid, selecNotValid), ...
                'CandidateSets:CandidateSets:InvalidArgument');
            
            selecNotValid = {true(1, 10), true(1, 10), true(1, 5)}';
            tc.assertError(@()CandidateSets (DataHandle (vertcat (candValid1, candValid2, candValid3)), lutValid, selecNotValid), ...
                'CandidateSets:CandidateSets:InvalidArgument');
            
            selecNotValid = {true(1, 10), true(1, 10), true(1, 5), true(1, 3)}';
            tc.assertError(@()CandidateSets (DataHandle (vertcat (candValid1, candValid2, candValid3)), lutValid, selecNotValid), ...
                'CandidateSets:CandidateSets:InvalidArgument');
            
            lutValid   = [1, 1, 2, NaN]';
            selecValid = {true(1, 10), true(1, 10), true(1, 5), NaN}';
            tc.assertInstanceOf (CandidateSets ( ...
                DataHandle (vertcat (candValid1, candValid2, candValid3)), ...
                lutValid, selecValid), 'CandidateSets');
            
            lutValid   = [1, 1, 2, NaN]';
            selecValid = {rand(1, 10) > 0.5, false(1, 10), true(1, 5), NaN}';
            tc.assertInstanceOf (CandidateSets ( ...
                DataHandle (vertcat (candValid1, candValid2, candValid3)), ...
                lutValid, selecValid), 'CandidateSets');
            
            % All structs should have the same fields
%             candValid1 = struct ('data', sparse(1, 10), 'id', {repmat({'a'}, [1, 10])}, 'num', 10, 'setId', 1);
%             candValid2 = struct ('id', {repmat({'a'}, [1, 5])}, 'data', ones(1, 5), 'num', 5);
%             candValid3 = struct ('data', 1, 'id', {'a'}, 'num', 1);
%             tc.assertError(@()CandidateSets (DataHandle ({candValid1, candValid2, candValid3}), lutValid), ...
%                 'CandidateSets:CandidateSets:InvalidArgument');
            
            % Some random scenarios
            for ii = 1:20000
                num1 = randi(10) - 1;
                num2 = randi(5) - 1;
                num3 = randi(3) - 1;
                candValid = struct ('data', [], 'id', {}, 'num', []);
                candValid(1) = struct ('data', rand (10, num1), 'id', {repmat({'a'}, [1, num1])}, 'num', num1);
                candValid(2) = struct ('id', {repmat({'b'}, [1, num2])}, 'data', sparse (10, num2), 'num', num2);
                candValid(3) = struct ('num', num3, 'id', {repmat({'c'}, [1, num3])}, 'data', rand (10, num3));

                lutValid = randi(4, randi(4) - 1, 1);
                makeNaNId = lutValid == 4;
                lutValid(makeNaNId) = NaN;

                selecValid = cell(size (lutValid));
                for i = 1:numel (lutValid)
                    if (isnan (lutValid(i)))
                        selecValid{i} = NaN;
                    else
                        selecValid{i} = rand (1, candValid(lutValid(i)).num) > 0.5;
                    end
                end % for

                if (any ([num1, num2, num3] == 0))
                    tc.assertError (@()CandidateSets ( ...
                        DataHandle (candValid), lutValid, selecValid), 'CandidateSets:CandidateSets:InvalidArgument');
                else
                    tc.assertInstanceOf (CandidateSets ( ...
                        DataHandle (candValid), lutValid, selecValid), 'CandidateSets');
                end % if
            end % for
            
        end % function       
        %% Test 3: Update of LUT and SELEC after "getSubset"
        function testGetSubset (tc)
            % create some dummy scenario
            d_Y = 10; 
            m = 20; 
            FP = rand(d_Y, m) > 0.5;
            INCHI = strsplit (sprintf ('inchi%d,', 1:m), ',');
            INCHI = INCHI(1:(end - 1));

            lut = [1, 1, 2, 2, 1, 3, 4, 3];
            
            cand = struct ('data', {}, 'id', {}, 'num', {});
            cand(1) = struct ('data', [FP(:, 1), FP(:, 7), FP(:, 3), FP(:, 19)],   ... 
                              'id', {{INCHI{1}, INCHI{7}, INCHI{3}, INCHI{19}}}, ...
                              'num', 4 );
            cand(2) = struct ('id', {{INCHI{2}, INCHI{5}, INCHI{20}, INCHI{18}, INCHI{12}}}, ...
                              'data', [FP(:, 2), FP(:, 5), FP(:, 20), FP(:, 18), FP(:, 12)], ... 
                              'num', 5 );                          
            cand(3) = struct ('num', 2 , 'data', [FP(:, 9), FP(:, 15)],   ... 
                              'id', {{INCHI{9}, INCHI{15}}});                          
            cand(4) = struct ('data', [FP(:, 1), FP(:, 11), FP(:, 13)],   ... 
                              'id', {{INCHI{1}, INCHI{11}, INCHI{13}}}, ...
                              'num', 3 );                          
                          
            % Add candidate structurs containing not informative fields
                   
            % lets take an empty subset
            Y_c = CandidateSets (DataHandle (cand), lut');
            tc.assertInstanceOf (Y_c.getSubset ([]), 'CandidateSets');
            Y_c_sub = Y_c.getSubset ([]);
            tc.assertEmpty (Y_c_sub.lut_);
            tc.assertEmpty (Y_c_sub.selec_);
            tc.assertTrue (Y_c_sub == CandidateSets (DataHandle (cand), []));
            tc.assertTrue (Y_c_sub == CandidateSets (DataHandle (cand)));
            clear Y_c Y_c_sub;
            
            % lets cause and out-of-range
            Y_c = CandidateSets (DataHandle (cand), lut');
            tc.assertInstanceOf (Y_c.getSubset ([]), 'CandidateSets');
            tc.assertError (@()Y_c.getSubset(-1),     'CandidateSets:getSubset:OutOfRange');
            tc.assertError (@()Y_c.getSubset(0),      'CandidateSets:getSubset:OutOfRange');
            tc.assertError (@()Y_c.getSubset(1:9),    'CandidateSets:getSubset:OutOfRange');
            tc.assertError (@()Y_c.getSubset(1:1000), 'CandidateSets:getSubset:OutOfRange');
            clear Y_c;
            
            Y_c = CandidateSets (DataHandle (cand), []);
            tc.assertInstanceOf (Y_c.getSubset ([]), 'CandidateSets');
            tc.assertError (@()Y_c.getSubset(1),      'CandidateSets:getSubset:OutOfRange');
            tc.assertError (@()Y_c.getSubset(-1),     'CandidateSets:getSubset:OutOfRange');
            tc.assertError (@()Y_c.getSubset(0),      'CandidateSets:getSubset:OutOfRange');
            tc.assertError (@()Y_c.getSubset(1:9),    'CandidateSets:getSubset:OutOfRange');
            tc.assertError (@()Y_c.getSubset(1:1000), 'CandidateSets:getSubset:OutOfRange');
            clear Y_c;
            
           
            % lets select all the elements and check preserved equalness
            Y_c = CandidateSets (DataHandle (cand), lut');
            tc.assertTrue (Y_c == Y_c.getSubset(1:8));
            clear Y_c;
            
            % lets select some more interesting subsets
            Y_c = CandidateSets (DataHandle (cand), lut');
            Y_c_sub = Y_c.getSubset (1:4);
            tc.assertTrue (all (size (Y_c.data_handle_.data_) == size (Y_c_sub.data_handle_.data_)));
            
            Y_c_sub_ref = CandidateSets (DataHandle (cand), [1, 1, 2, 2]');
            tc.assertTrue (Y_c_sub == Y_c_sub_ref);
            clear Y_c_sub Y_c_sub_ref;
            
            Y_c_sub = Y_c.getSubset ([5, 1, 8]);
            Y_c_sub_ref = CandidateSets (DataHandle (cand), lut([5, 1, 8])');
            tc.assertTrue (Y_c_sub == Y_c_sub_ref);
            clear Y_c_sub Y_c_sub_ref;
            
            % Get subset of subset !!!
            Y_c = CandidateSets (DataHandle (cand), lut');
            Y_c_sub = Y_c.getSubset ([5, 1, 8]);
            Y_c_sub_sub = Y_c_sub.getSubset ([1, 3]);
            tc.assertTrue (Y_c_sub_sub == CandidateSets ( ...
                DataHandle (cand), lut([5, 8])'));
            
            selec_fake = cell (2, 1); 
            selec_fake{1} = logical ([1, 1, 0, 1]);
            selec_fake{2} = logical ([1, 1]);
            tc.assertFalse (Y_c_sub_sub == CandidateSets ( ...
                DataHandle (cand), lut([5, 8])', selec_fake));      
        end % function 
        %% Test 4: isEqual
        function testIsEqual (tc)
            % create some dummy scenario
            d_Y = 10; 
            m = 20; 
            FP = rand(d_Y, m) > 0.5;
            INCHI = strsplit (sprintf ('inchi%d,', 1:m), ',');
            INCHI = INCHI(1:(end - 1));

            lut = [1, 1, 2, 2, 1, 3, 4, 3];
            
            cand = struct([]);
            cand(1).data = [FP(:, 1), FP(:, 7), FP(:, 3), FP(:, 19)]; 
            cand(1).id   = {INCHI{1}, INCHI{7}, INCHI{3}, INCHI{19}};
            cand(1).num  = 4;
            
            cand(2) = struct ('data', [FP(:, 2), FP(:, 5), FP(:, 20), FP(:, 18), FP(:, 12)],   ... 
                              'id', {{INCHI{2}, INCHI{5}, INCHI{20}, INCHI{18}, INCHI{12}}}, ...
                              'num', 5 );                          
            cand(3) = struct ('data', [FP(:, 9), FP(:, 15)],   ... 
                              'id', {{INCHI{9}, INCHI{15}}}, ...
                              'num', 2 );                          
            cand(4) = struct ('data', [FP(:, 1), FP(:, 11), FP(:, 13)],   ... 
                              'id', {{INCHI{1}, INCHI{11}, INCHI{13}}}, ...
                              'num', 3 );                         
            
            % Lets test whether two empty elements are equal                          
            Y_c = CandidateSets (DataHandle (cand(1)), []);
            tc.assertTrue (Y_c == Y_c);
            clear Y_c:
                          
            % Lets test whether two times the same element are equal
            Y_c = CandidateSets (DataHandle (cand), lut');
            tc.assertTrue (Y_c == Y_c);
            clear Y_c:
            
            % Lets test whether two a candidate set is equal with the
            % candidate set which is the full-subset of it
            Y_c = CandidateSets (DataHandle (cand), lut');
            tc.assertTrue (Y_c == Y_c.getSubset(1:8));
            clear Y_c;
            
            % Lets test some not equal things
            Y_c   = CandidateSets (DataHandle (cand), lut');
            Y_c_2 = CandidateSets (DataHandle (cand), [1, 2]');
            tc.assertFalse (Y_c == Y_c_2);
            clear Y_c_2;
            
            tc.assertFalse (Y_c == 5);
            
            cand2 = cand; 
            cand2(1).data(3, 2) = ~cand2(1).data(3, 2);
            Y_c_2 = CandidateSets (DataHandle (cand2), lut');
            tc.assertFalse (Y_c == Y_c_2);
            clear cand2 Y_c_2;
            
            cand2 = cand; 
            cand2(1).id{1} = 'FAQ';           
            Y_c_2 = CandidateSets (DataHandle (cand2), lut');
            tc.assertFalse (Y_c == Y_c_2);
            clear cand2 Y_c_2;
            
            % This test requires to set the object's properties public
%             cand2 = cand; 
% %             cand2{2}.num = 1000;           
%             Y_c_2 = CandidateSets (DataHandle (cand2), lut');
%             Y_c_2.data_handle_.data_{2}.num = 1000;
%             tc.assertFalse (Y_c == Y_c_2);
%             clear cand2 Y_c_2;
%             clear Y_c;
%             
            % lut = [1, 1, 2, 2, 1, 3, 4, 3];
            Y_c   = CandidateSets (DataHandle (cand), lut');
            selec = cell (4, 1);
            selec{1} = logical([1, 1, 1, 1]);
            selec{2} = logical([1, 1, 1, 1]);
            selec{3} = logical([1, 1, 1, 1, 1]);
            selec{4} = logical([1, 1, 1, 1, 1]);
            selec{5} = logical([1, 1, 1, 1]);
            selec{6} = logical([1, 1]);
            selec{7} = logical([1, 1, 1]);
            selec{8} = logical([1, 1]);
            Y_c_2 = CandidateSets (DataHandle (cand), lut', selec);
            tc.assertTrue (Y_c == Y_c_2);
            
            selec{1} = logical([1, 1, 0, 1]);
            Y_c_2 = CandidateSets (DataHandle (cand), lut', selec);
            tc.assertFalse (Y_c == Y_c_2);    
            
            % CandidateSets are equal if the candiate set data structures
            % contain the same fieldnames. However, the cotent of those 
            % fields is only compared if this field is one of the standards
            % 'data', 'id' and 'num'
            cand2 = struct ([]);
            
            cand2(1).data  = [FP(:, 1), FP(:, 7), FP(:, 3), FP(:, 19)];
            cand2(1).id    = {INCHI{1}, INCHI{7}, INCHI{3}, INCHI{19}};
            cand2(1).num   = 4;
            cand2(1).setId = 1;
            
            cand2(2) = struct ('data', [FP(:, 2), FP(:, 5), FP(:, 20), FP(:, 18), FP(:, 12)],... 
                              'id', {{INCHI{2}, INCHI{5}, INCHI{20}, INCHI{18}, INCHI{12}}}, ...
                              'num', 5, 'setId', 2);                          
            cand2(3) = struct ('data', [FP(:, 9), FP(:, 15)],   ... 
                              'id', {{INCHI{9}, INCHI{15}}}, ...
                              'setId', 3, 'num', 2);                          
            cand2(4) = struct ('num', 3, 'data', [FP(:, 1), FP(:, 11), FP(:, 13)],   ... 
                              'id', {{INCHI{1}, INCHI{11}, INCHI{13}}}, ...
                              'setId', 4);  
                          
            % The class of the "main" field should be the same
            cand_new_class = cand2;
            cand_new_class(1).data = sparse (10, 4);
            Y_c = CandidateSets (DataHandle (cand2), lut', selec);
            Y_c_2 = CandidateSets (DataHandle (cand_new_class), lut', selec);
            tc.assertFalse (Y_c == Y_c_2);     
            
            % fieldnames are not the same
            Y_c   = CandidateSets (DataHandle (cand), lut', selec);
            Y_c_2 = CandidateSets (DataHandle (cand2), lut', selec);
            tc.assertFalse (Y_c == Y_c_2);     
            
            % If the fieldnames are the same, only the three main fields
            % are important.
            Y_c   = CandidateSets (DataHandle (cand2), lut', selec);
            cand2(4).setId = 99;
            Y_c_2 = CandidateSets (DataHandle (cand2), lut', selec);
            tc.assertTrue (Y_c == Y_c_2);                 
        end % function
        
        %% Test 5: getCandidateSet
        function testGetCandidateSet (tc)
            % create some dummy scenario
            d_Y = 10; 
            m = 20; 
            FP = rand(d_Y, m) > 0.5;
            INCHI = strsplit (sprintf ('inchi%d,', 1:m), ',');
            INCHI = INCHI(1:(end - 1));

            lut = [1, 1, 2, 2, 1, 3, 4, 3];
            
            cand = struct ('data', {}, 'id', {}, 'num', {}, 'setId', {});
            cand(1) = struct ('data', sparse([FP(:, 1), FP(:, 7), FP(:, 3), FP(:, 19)]),   ... 
                              'id', {{INCHI{1}, INCHI{7}, INCHI{3}, INCHI{19}}}, ...
                              'num', 4, 'setId', 1);
            cand(2) = struct ('num', 5, 'setId', 2, ...
                              'data', [FP(:, 2), FP(:, 5), FP(:, 20), FP(:, 18), FP(:, 12)],   ... 
                              'id', {{INCHI{2}, INCHI{5}, INCHI{20}, INCHI{18}, INCHI{12}}});                          
            cand(3) = struct ('data', [FP(:, 9), FP(:, 15)],   ... 
                              'id', {{INCHI{9}, INCHI{15}}}, ...
                              'num', 2, 'setId', 3);                          
            cand(4) = struct ('id', {{INCHI{1}, INCHI{11}, INCHI{13}}}, ...
                              'data', [FP(:, 1), FP(:, 11), FP(:, 13)], ...
                              'num', 3, 'setId', 4);    
                          
            % Let the function fail because of malformed input
            Y_c = CandidateSets (DataHandle (cand), lut');
            tc.assertError (@()Y_c.getCandidateSet([1, 2]), 'CandidateSets:getCandidateSet:InvalidArgument');
            tc.assertEmpty (Y_c.getCandidateSet([]));
            clear Y_c;
            
            % lets get some elements 
            Y_c = CandidateSets (DataHandle (cand), lut');
            for i = 1:numel (lut)
                candSet_i = Y_c.getCandidateSet (i);
                tc.assertEqual (candSet_i.setId, lut(i));
            end % if   
            clear Y_c candSet_i;
            
            % lets check whether a nan candset is removed if needed
            lut2 = lut; 
            lut2(2) = NaN;
            Y_c = CandidateSets (DataHandle (cand), lut2');
            candSet_i = Y_c.getCandidateSet (2);
            tc.assertTrue (isnan (candSet_i));         
            
            % lets get some elements using selection
            Y_c = CandidateSets (DataHandle (cand), lut');
            for i = 1:numel (lut)
                candSet_i = Y_c.getCandidateSet (i, 1);
                tc.assertEqual (candSet_i, cand(lut(i)));
            end % if   
            clear Y_c candSet_i;
                   
            selec = cell (8, 1);
            selec{1} = logical([0, 1, 1, 1]);
            selec{2} = logical([1, 1, 1, 1]);
            selec{3} = logical([1, 1, 0, 1, 1]);
            selec{4} = logical([1, 1, 1, 1, 1]);
            selec{5} = logical([1, 1, 1, 1]);
            selec{6} = logical([1, 0]);
            selec{7} = logical([1, 1, 1]);
            selec{8} = logical([1, 1]);
            Y_c = CandidateSets (DataHandle (cand), lut', selec);
            for i = 1:numel (lut)
                candSet_i = Y_c.getCandidateSet (i, 1);
                
                tc.assertEqual (candSet_i.num, sum(selec{i}));
                tc.assertTrue (all (size (candSet_i.id) == size (cand(lut(i)).id(selec{i}))));
                tc.assertTrue (all (size (candSet_i.data) == size (cand(lut(i)).data(:, selec{i}))));
                tc.assertTrue (candSet_i.setId == cand(lut(i)).setId);
                
                candTmp = struct ('id', {cand(lut(i)).id(selec{i})}, ...
                    'data', cand(lut(i)).data(:, selec{i}), ...
                    'num', sum (selec{i}), ...
                    'setId', cand(lut(i)).setId);
                
                tc.assertEqual (candSet_i, candTmp);
            end % if   
            clear Y_c candSet_i;
            
            selec = cell (8, 1);
            selec{1} = logical([0, 0, 0, 0]);
            selec{2} = logical([1, 1, 1, 1]);
            selec{3} = logical([1, 1, 0, 1, 1]);
            selec{4} = logical([1, 1, 1, 1, 1]);
            selec{5} = logical([1, 1, 1, 1]);
            selec{6} = logical([0, 0]);
            selec{7} = logical([1, 1, 1]);
            selec{8} = logical([1, 1]);
            Y_c = CandidateSets (DataHandle (cand), lut', selec);
            for i = 1:numel (lut)
                candSet_i = Y_c.getCandidateSet (i, 1);
                
                tc.assertEqual (candSet_i.num, sum(selec{i}));
                tc.assertTrue (all (size (candSet_i.id) == size (cand(lut(i)).id(selec{i}))));
                tc.assertTrue (all (size (candSet_i.data) == size (cand(lut(i)).data(:, selec{i}))));
                tc.assertTrue (candSet_i.setId == cand(lut(i)).setId);
                
                % Even if 'setId' does not belong to the main-fields it
                % needs to stay after the selection
                candTmp = struct ('id', {cand(lut(i)).id(selec{i})}, ...
                    'data', cand(lut(i)).data(:, selec{i}), ...
                    'num', sum (selec{i}), ...
                    'setId', cand(lut(i)).setId);
                tc.assertEqual (candSet_i, candTmp);
                
                candTmp = struct ('id', {cand(lut(i)).id(selec{i})}, ...
                    'data', cand(lut(i)).data(:, selec{i}), ...
                    'num', sum (selec{i}));
                tc.assertNotEqual (candSet_i, candTmp);
            end % if   
            clear Y_c candSet_i;
            
            % Add some randomized tests here
            newFields = {'bla', 'blub', 'huhu'};
            illigalField = {'nixhier'};
            
            for run_i = 1:20000
                % Create random candidate set
                num1 = randi(10);
                num2 = randi(5);
                num3 = randi(3);
                candValid = struct ('data', {}, 'id', {}, 'num', {}, 'setId', {});
                candValid(1) = struct ('data', rand (20, num1), 'id', {repmat({'a'}, [1, num1])}, 'setId', rand, 'num', num1);
                candValid(2) = struct ('id', {repmat({'b'}, [1, num2])}, 'data', sparse (20, num2), 'num', num2, 'setId', rand);
                candValid(3) = struct ('num', num3, 'id', {repmat({'c'}, [1, num3])}, 'data', rand (20, num3), 'setId', rand);
                
                % Create a random LUT
                lut = randi(4, randi(20) - 1, 1);
                makeNaNId = lut == 4;
                lut(makeNaNId) = NaN;

                % Create a random SELECTION
                selec = cell(size (lut));
                for i = 1:numel (lut)
                    if (isnan (lut(i)))
                        selec{i} = NaN;
                    else
                        selec{i} = rand (1, candValid(lut(i)).num) > 0.5;
                    end
                end % for
                
                % Add some random fields to see, whether they are preserved
                % without any change 
                n_newFields  = randi (numel (newFields));
                id_newFields = randperm (numel (newFields));
                id_newFields = id_newFields(1:n_newFields);
                
                for i = id_newFields ; for j = 1:numel (candValid)
                    candValid(j).(newFields{i}) = rand(9, 1);                     
                end ; end % for
            
                % Build up the candidate set object
                Y_c = CandidateSets (DataHandle (candValid), lut, selec);
                
                % Should we take only the selected candidates?
                doSelection = rand > 0.5;
                
                % Choose a random field to extract
                availableFields = horzcat ({'data', 'id', 'num', 'setId'}, newFields(id_newFields), illigalField);
                id_field = randi(numel (availableFields));
                fieldSelection = availableFields{id_field};
                
                % Get the candidate set corresponding to a random id
                if (isempty (lut))
                    candSet_id = Y_c.getCandidateSet ([], doSelection);
                    field_id = Y_c.getCandidateSet ([], doSelection, fieldSelection);
                    
                    tc.assertEmpty (candSet_id);
                    tc.assertEmpty (field_id);
                    continue;
                else
                    id = randi (length (lut));
                    candSet_id = Y_c.getCandidateSet (id, doSelection);
                    
                    % If there is no candidate set for this id (NaN)
                    if (isnan (lut(id)))                                            
                        field_id = Y_c.getCandidateSet (id, doSelection, fieldSelection);
                        
                        tc.assertTrue (isnan (candSet_id));
                        tc.assertTrue (isnan (field_id));

                        continue;
                    else
                        if (strcmp (fieldSelection, illigalField))
                            tc.assertError (@()Y_c.getCandidateSet (id, doSelection, fieldSelection), ...
                                'CandidateSet:getCandidateSet:InvalidArgument');

                            continue;
                        else
                            field_id = Y_c.getCandidateSet (id, doSelection, fieldSelection);
                        end % if                         
                    end % if
                end % if 


                % If a selection has been done and this selection was
                % empty:
                if (all (selec{id} == false) && doSelection)
                    tc.assertEmpty (candSet_id.id);
                    tc.assertEmpty (candSet_id.data);
                    tc.assertTrue (candSet_id.num == 0);
                    
                    switch (fieldSelection)
                        case {'data', 'id'}
                            tc.assertEmpty (field_id);
                        case 'num'
                            tc.assertTrue (field_id == 0);
                        otherwise
                            tc.assertTrue (all (size (field_id) == size (candValid(lut(id)).(fieldSelection))));
                            tc.assertTrue (strcmp (class (field_id), class (candValid(lut(id)).(fieldSelection))));
                    end % switch
                    
                    continue;
                end % if
                
                if (doSelection)
                    % Check size
                    tc.assertTrue (all (size (candSet_id.id) == size (candValid(lut(id)).id(selec{id}))));
                    tc.assertTrue (all (size (candSet_id.data) == size (candValid(lut(id)).data(:, selec{id}))));
                    
                    % Check content
                    tc.assertTrue (candSet_id.setId == candValid(lut(id)).setId);
                    tc.assertTrue (candSet_id.num == sum(selec{id}));
                    tc.assertTrue (all (all (candSet_id.data == candValid(lut(id)).data(:, selec{id}))));
                    tc.assertTrue (all (strcmp (candSet_id.id, candValid(lut(id)).id(selec{id}))));
                    
                    % Check only the selected field
                    switch (fieldSelection)
                        case 'data'
                            tc.assertTrue (all (size (field_id) == size (candValid(lut(id)).data(:, selec{id}))));
                            tc.assertTrue (all (all (field_id == candValid(lut(id)).data(:, selec{id}))));
                        case 'id'
                            tc.assertTrue (all (size (field_id) == size (candValid(lut(id)).id(selec{id}))));
                            tc.assertTrue (all (strcmp (field_id, candValid(lut(id)).id(selec{id}))));
                        case 'num'
                            tc.assertTrue (field_id == sum(selec{id}));
                        otherwise 
                            % Check only the selected field
                            tc.assertTrue (all (size (field_id) == size (candValid(lut(id)).(fieldSelection))));
                            tc.assertTrue (strcmp (class (field_id), class (candValid(lut(id)).(fieldSelection))));
                    end % switch       
                else
                    % Check size
                    tc.assertTrue (all (size (candSet_id.id) == size (candValid(lut(id)).id)));
                    tc.assertTrue (all (size (candSet_id.data) == size (candValid(lut(id)).data)));
                    
                    % Check content 
                    tc.assertEqual (candSet_id.num, candValid(lut(id)).num);
                    tc.assertTrue (candSet_id.setId == candValid(lut(id)).setId);
                    tc.assertTrue (all (all (candSet_id.data == candValid(lut(id)).data)));
                    tc.assertTrue (all (strcmp (candSet_id.id, candValid(lut(id)).id)));
                    
                    % Check only the selected field
                    switch (fieldSelection)
                        case 'data'
                            tc.assertTrue (all( all (field_id == candValid(lut(id)).data)));
                        case 'id'
                            tc.assertTrue (all (strcmp (field_id, candValid(lut(id)).id)));
                        case 'num'
                            tc.assertTrue (field_id == candValid(lut(id)).num);
                        otherwise 
                            ;
                    end % switch  
                    
                    % Check only the selected field
                    tc.assertTrue (all (size (field_id) == size (candValid(lut(id)).(fieldSelection))));
                    tc.assertTrue (strcmp (class (field_id), class (candValid(lut(id)).(fieldSelection))));
                end % if
                
                % Check additional fields' content
                tc.assertTrue (all (ismember (newFields(id_newFields), fieldnames (candSet_id))));
                for i = id_newFields
                    tc.assertEqual (candSet_id.(newFields{i}), candValid(lut(id)).(newFields{i}));
                end % for
            end % for        
        end % function       
    end % methods
end % class
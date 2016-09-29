function ranks = getRanksBasedOnScores (Y_C, inchis, scores, eval_set)
%% GETRANKSBASEDONSCORES calculates the pr-image rank for each example 
%    This function used the prediced scores for each candidate given a
%    certain test-example. The canidates are sorted according to their
%    score. Subsequently the test-example is searched in the set of sorted
%    candidates. The position at which it is found is it's rank.

    ranks = NaN (Y_C.getNumberOfExamples(), 1);

    for j = 1:Y_C.getNumberOfExamples()
        % In the case of the metabolite identification we should only
        % evaluate the examples from the eval-set. This ensures
        % comparable results with the recent publication of Celine
        % [Brouard2016].
        if (~ eval_set(j)) ; continue; end % if

        % There might be scores which are NaN, as not each test example
        % has a candidate set. 
        assert (~ any (isnan (scores{j})), 'We only evaluate using the eval-set.');
        assert (~ isnan (Y_C.getCandidateSet (j, 0, 'num')), ...
            'We only evaluate using the eval-set.');

        [~ , IX] = sort (scores{j}, 'descend');

        % Get the inchis of all the candidate in the set for test
        % example j
        inchis_c = Y_C.getCandidateSet (j, 0, 'id');
        ind = find (strcmp (inchis_c(IX), inchis(j)));
        assert (~ isempty (ind), ...
            'The inchi of the example itself is not in the candidate set. We only evaluate using the eval-set.');

        ranks(j) = find (strcmp (inchis_c(IX), inchis{j}));
    end % for
end % function
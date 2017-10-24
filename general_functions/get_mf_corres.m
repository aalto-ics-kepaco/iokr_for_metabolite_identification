function mf_corres = get_mf_corres (mfs, cand)
    assert (iscell (mfs));

    mfs_cand = arrayfun (@(x) x.set_desc, cand, 'UniformOutput', false);
    [~, mf_corres] = ismember (mfs, mfs_cand);
    
    if ~ (all (mf_corres > 0))
        warning ('get_mf_corres:InvalidInput', 'Only %d molecular formulas out of %d are in the candidate set.', ...
            sum (mf_corres > 0), length (mfs));
    end % if
end % function



function mf_corres = get_mf_corres (mfs, cand)
    mfs_cand = arrayfun (@(x) x.set_desc, cand, 'UniformOutput', false);
    [~, mf_corres] = ismember (mfs, mfs_cand);
    assert (all (mf_corres > 0));
end % function
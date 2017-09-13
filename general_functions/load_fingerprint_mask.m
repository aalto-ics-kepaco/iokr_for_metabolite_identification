function fp_mask_bin = load_fingerprint_mask (fp_mask_fn)

    fp_mask_str = strsplit (fileread (fp_mask_fn), '\t');
    % x:      property is used
    % t:      property is always true
    % f:      property is always false
    % NUMBER: property is identical with the one on bit-position NUMBER
    fp_mask_bin = cellfun (@(c) all (c == 'x'), fp_mask_str, ...
        'UniformOutput', true);
    
end % function

function prepare_data_for_iokr (wdir)
    % Load data
    inchis  = readtable ([wdir, '/ms/compound_list'], ...
        'ReadVariableNames', false, 'Delimiter', '\t');
    inchis.Properties.VariableNames = {'mol_id', 'inchi_key_1', 'inchi'};
    mfs     = readtable ( [wdir, '/ms/ml_formula_per_compound'], ...
        'ReadVariableNames', false, 'Delimiter', ' ');
    mfs.Properties.VariableNames = {'mol_id', 'molecular_formula'};
    dt_inchi_mf = join (inchis, mfs, 'Keys', 'mol_id');

    load ( [wdir, '/fingerprints/fingerprints.mat'], 'dt_fp');
    
    % Restrict inchis and molecular formulas to those present in the kernel 
    % matrix.
    %dt_inchi_mf = dt_inchi_mf( ...
    %    cellfun (@(c) ismember (c, mol_ids), dt_inchi_mf.mol_id), :);
    %assert (all (strcmp (mol_ids, dt_inchi_mf.mol_id)));

    % Build table containing the:
    % mol_id:            unique identifier for each training spectra
    % inchi_key_1:       unique identifier for each training compound
    % fp_full:           full fingerprint for each training compound
    % fp_masked:         masked fingerprint for each training compound
    % molecular_formula: molecular formula for each training compound
    dt_inchi_mf_fp = join (dt_inchi_mf, dt_fp, 'Keys', 'mol_id');

    save ([wdir, '/compounds/compound_info.mat'], 'dt_inchi_mf_fp', '-v7.3');
    % Check the candidate set and add missing true compounds
    load ([wdir, '/candidates/cand.mat']);
    
    % Is there a candidate set for each training compound?
    mfs_cand = arrayfun (@(x) x.set_desc, cand, 'UniformOutput', false);
    assert (all (ismember (dt_inchi_mf_fp.molecular_formula, mfs_cand)), ...
        'Missing candidate set.');

    non_matching_fps = 0;
    
    for i_cmp = 1:size (dt_inchi_mf_fp, 1)
        fprintf ('Process each compounds %d/%d\n', i_cmp, size (dt_inchi_mf_fp, 1));

        % Is each compound in the candidate?
        inchi_key_cmp = char (dt_inchi_mf_fp{i_cmp, 'inchi_key_1'});
        mf_cmp        = char (dt_inchi_mf_fp{i_cmp, 'molecular_formula'});
        fp_cmp        = dt_inchi_mf_fp{i_cmp, 'fp_full'};
        assert (issparse (fp_cmp));

        i_cand = find (strcmp (mfs_cand, mf_cmp)); 
        assert (length (i_cand) == 1);

        % If the current inchi_key is not in the candidate set, than put the it
        % into the set and use the fingerprint from the training.

        [~, locb] = ismember (inchi_key_cmp, cand(i_cand).id);
        if locb == 0
            warning ('True compound is not in the candidate set.');

            cand(i_cand).num  = cand(i_cand).num + 1;
            cand(i_cand).id   = [cand(i_cand).id; inchi_key_cmp];
            cand(i_cand).data = [cand(i_cand).data, fp_cmp'];
        else 
            if (~ all (cand(i_cand).data(:, locb) == fp_cmp'))
                warning ('Fingerprints do not match.');

                disp (cand(i_cand).data(:, locb) ~= fp_cmp')
                
                non_matching_fps = non_matching_fps + 1;
                
                cand(i_cand).data(:, locb) = fp_cmp';
            end % if
        end % if 
    end % for
    
    for i_cand = 1:length (cand)
        fprintf ('%d/%d\n', i_cand, length (cand));

        cand(i_cand).id = cand(i_cand).id';
    end % for
    
    disp (non_matching_fps)

    save ([wdir, '/candidates/cand_prepared.mat'], 'cand', '-v7.3');
end % function
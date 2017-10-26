function prepare_compound_info_for_iokr (wdir)
%% PREPARE_COMPOUND_INFO_FOR_IORK
% Build table containing the:
% mol_id:            unique identifier for each training spectra
% inchi_key_1:       unique identifier for each training compound
% molecular_formula: molecular formula for each training compound
% fp_full:           full fingerprint for each training compound
% fp_masked:         masked fingerprint for each training compound

    if (~ exist ([wdir, '/fingerprints/fingerprints.mat'], 'file'))
        prepare_fingerprints_for_iokr (wdir);
    end % if

    % Load compound list
    cmps  = readtable ([wdir, '/spectra/compound_list'], ...
        'ReadVariableNames', false, 'Delimiter', '\t');
    cmps.Properties.VariableNames = { ...
        'mol_id', ...
        'inchi_key_1', ...
        'inchi', ...
        'molecular_formula'};
    
    % Load fingerprints
    load ([wdir, '/fingerprints/fingerprints.mat'], 'dt_fp');
 
    dt_inchi_mf_fp = join (cmps, dt_fp, 'Keys', 'mol_id');

    save ([wdir, '/spectra/compound_info.mat'], 'dt_inchi_mf_fp', '-v7.3');
end % function
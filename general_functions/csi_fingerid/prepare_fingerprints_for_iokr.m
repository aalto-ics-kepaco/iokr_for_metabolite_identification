function prepare_fingerprints_for_iokr (wdir)
%% PREPARE_FINGERPRINTS_FOR_IOKR read fingerprint files from csi-fingerid
%   Csi-fingerid stores the fingerprints for each compound in a separate
%   file. This script reads all files in the provided directory. The
%   fingerprints are stored in a DataTable with 3 columns:
%       'mol_id': String based on the filename
%       'fp_full': Sparse vector of all fingerprints
%       'fp_masked': Sparse vector of the masked fingerprints

    fp_files = dir ([wdir, '/fingerprints/*.fpt']);

    fp_mask_bin = load_fingerprint_mask ([wdir, '/fingerprints/fingerprints.mask']);

    dt_fp = cell2table (cell (0, 3), ...
        'VariableNames', {'mol_id', 'fp_full', 'fp_masked'});

    for fp_file = fp_files'
        % Open file
        file_id = fopen ([wdir, '/fingerprints/', fp_file.name], 'r');

        % Read fingerprint string
        fp_str = textscan (file_id, '%s');
        fp_str = char (fp_str{1, 1});

        % Convert string into vector
        fp        = sparse (fp_str - '0');
        fp_masked = fp (fp_mask_bin);

        % Add fingerprint to table
        mol_id = basename (fp_file.name);
        dt_fp = [dt_fp; {mol_id, fp, fp_masked}];

        % Close file
        fclose (file_id);
    end % for

    save ([wdir, '/fingerprints/fingerprints.mat'], 'dt_fp', '-v7.3');
end % function 
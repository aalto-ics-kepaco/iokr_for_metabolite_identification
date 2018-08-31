function run_IOKR_independent_set(input_dir_training, input_dir_indep, output_dir, iokr_param)
%% RUN_IOKR_INDEPENDENT_SET Scoring script for independent set of spectra
% Requires preprocessed data: path/to/independent_set
%   - kernels/      : kernel files between independent MSMS and training
%                     data
%   - candidates/   : csv-files containing candidates, fingerprints for all
%                     candidates

    %--------------------------------------------------------------
    % Load and prepare data
    %--------------------------------------------------------------

    % Load the list of training compounds. This list can be created using
    % the csi-fingerid CLI: 'fingerID list-compounds > compounds'
    cmps_train = loadCompoundList(input_dir_training, 'csi_fingerid');
    
    % Load the list of independent compounds. This list can be created
    % using the csi-fingerid CLI: 'fingerID list-independent-compounds --indep-set=INDEPSET > independent/INDEPSET/compounds'
    % NOTE: If the compounds are not known, i.e. no InChI or SMILES is
    %       available in the ms-file. Than the table contains only the
    %       'spec_id', e.g. the 'inchikey2D' column would be filled with
    %       '__NOT__KNOWN__'.
    cmps_indep = loadCompoundList(input_dir_indep, 'csi_fingerid');
    
    
end % function


% <<<<<<< HEAD
% %% Load list of molecular formulas
% mls_fn  = '/home/bach/Documents/studies/doctoral/data/csi_fingerid/fingerid-110/spectra/ml_formula_unique_all.txt';
% mls_str = table2array (readtable (mls_fn, 'ReadVariableNames', false));
% 
% n_mls = length (mls_str);
% fprintf ('Number of (unique) molecular formula: %d\n', n_mls);
% 
% 
% %% Load the fingerprint mask
% fp_mask_fn  = '/home/bach/Documents/studies/doctoral/data/csi_fingerid/fingerid-110/fingerprints/fingerprints.mask';
% fp_mask_bin = load_fingerprint_mask (fp_mask_fn);
% 
% % Determine the length of the stored 
% n_fps     = length (fp_mask_bin);
% n_fps_set = sum (fp_mask_bin);
% 
% fprintf ('Fingerprints: %d of %d bits are used.\n', ...
%     n_fps_set, n_fps);
% 
% %% Connect to the pubchem database
% % Load the driver
% % javaaddpath ('/usr/share/java/postgresql.jar')
% 
% % Database login, location, etc ...
% % datasource = 'pubchem';
% % username   = 'fingerid';
% % password   = 'tV9QRQHn2THjq5HR';
% % driver     = 'org.postgresql.Driver';
% % url        = 'jdbc:postgresql://pubchem.bioinf.uni-jena.de:5432/pubchem';
% 
% datasource = 'pubchem';
% username   = 'eric';
% password   = '47ex3z9tA';
% driver     = 'org.postgresql.Driver';
% url        = 'jdbc:postgresql://localhost:5433/pubchem';
% 
% % Establish connection
% conn = database (datasource, username, password, driver, url, ...
%     'ReadOnly','on');
% 
% if ~ isopen (conn)
%     error ('Could not establish db-connection: %s', ...
%         conn.Message);
% end % if
% 
% %% Download the canidate sets from the pubchem database
% % Create candidate set structure
% % id:     Identifier for each molecule
% % data:   Fingerprint for each molecule
% % dim:    Dimension of the fingerprints
% % num:    Number of molecules in the candidate set
% % set_id: Molecular formula used to create the candidate set
% cand_FN = '/home/bach/Documents/studies/doctoral/data/csi_fingerid/fingerid-110/fingerid-110/candidates/cand.mat';
% 
% % cand = repmat (struct ('id', NaN, 'data', NaN, ...
% %     'dim', NaN, 'num', NaN, 'set_desc', NaN), ...
% %     n_mls, 1);
% % save (cand_FN, 'cand', '-v7.3');
% % clear cand;
% 
% mfile = matfile (cand_FN, 'Writable', true);
% 
% i_ml_failed = [];
% 
% for i_ml = 206:n_mls
%     tic;
%     fprintf ('Process molecular formula: %s (%d/%d)\n', ...
%         mls_str{i_ml}, i_ml, n_mls);
%     
%     selec_str = sprintf ( ...
%         'SELECT inchi_key_1, fingerprint FROM fingerprints WHERE formula = ''%s'' AND fp_id = 1;', ...
%         mls_str{i_ml});
%     % Query DB
%     %tic;
%     cur     = exec (conn, selec_str);
%     data_db = fetch (cur);
%     %toc;
%     
%     % Parse information of the current candidate set
%     n_cand = size (data_db.Data, 1);
%     fprintf ('Number of candidates: %d\n', n_cand);
% =======
function conn = connect_to_pubchem_db (userdata)
    %% Connect to the pubchem database
% >>>>>>> e08abdc64dfab16861b3c3eaa22ac36952bdd6db
    
    % Load the driver
    % javaaddpath ('/u/33/bache1/unix/lib/postgresql-jdbc3-9.2.jar')
    javaaddpath ('/usr/share/java/postgresql.jar')
    
    % Database login, location, etc ...
    % datasource = 'pubchem';
    % username   = 'fingerid';
    % password   = 'tV9QRQHn2THjq5HR';
    % driver     = 'org.postgresql.Driver';
    % url        = 'jdbc:postgresql://pubchem.bioinf.uni-jena.de:5432/pubchem';

    datasource = 'pubchem';
    username   = userdata.username;
    password   = userdata.password;
    driver     = 'org.postgresql.Driver';
    url        = 'jdbc:postgresql://localhost:5433/pubchem';

    % Establish connection
    conn = database (datasource, username, password, driver, url, ...
        'ReadOnly','on');

    if ~ isopen (conn)
        error ('Could not establish db-connection: %s', ...
            conn.Message);
    end % if
end % function 
function fps = loadFingerprints(fps_dir, spec_ids, verbose)
%% LOADFINGERPRINTS Load fingerprints from *.fpt's or *.mat file
    if nargin < 2
        spec_ids = {};
    end % if
    if nargin < 3
        verbose = true;
    end % if

    fps_mat_fn_basename = "fps.mat";
    fps_mat_fn = strcat(fps_dir, "/", fps_mat_fn_basename);
    if exist(fps_mat_fn, 'file')  
        if verbose 
            fprintf('Load "fps.mat" file.\n')
        end % if
        
        fps = load(fps_mat_fn);
        fps = fps.fps;
    else
        n_spec = length(spec_ids);
        
        if verbose
            fprintf('Load separated *.fpt files: # = %d. Create "fps.mat" subsequently.\n', n_spec);
        end % if 
        
        % Load the fingerprint mask and determine fingerprint dimension
        fps_mask = load_fingerprint_mask(strcat(fps_dir, "/fingerprints.mask"));
        n_fps = length(fps_mask);
        
        fps = struct('fps', sparse(n_fps, n_spec), 'mask', fps_mask);
        
        for ii = 1:n_spec
            if verbose &&  ~ mod(ii, 500)
                fprintf('Process spectra id %d.\n', ii);
            end % if
            
            fps_spec_fn = strcat(fps_dir, "/", spec_ids{ii}, ".fpt");
            fid = fopen(fps_spec_fn, 'r');
            if fid == -1
                error('loadFingerprints:NoSuchFileOrDirectory', 'Cannot open fps-file "%s".', ...
                    fps_spec_fn);
            end % if
            
            % parse fingerprint vector 
            fps.fps((fgetl(fid) - '0') == 1, ii) = true;
            
            fclose(fid);
        end % if
        
        % save the fingerprint structure fur future use
        save(fps_mat_fn, 'fps', '-v7.3');
    end % if
end % function
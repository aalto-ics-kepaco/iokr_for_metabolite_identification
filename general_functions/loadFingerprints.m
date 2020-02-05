function [fps, mask] = loadFingerprints(raw_fps_dir, mat_fps_dir, spec_ids, verbose)
%% LOADFINGERPRINTS Load fingerprints from *.fpt's or *.mat file
    if nargin < 2
        mat_fps_dir = raw_fps_dir;
    end % if 
    if nargin < 3
        spec_ids = {};
    end % if
    if nargin < 4
        verbose = true;
    end % if

    % Calculate a hash value for the specified spec_ids so that the
    % fps_HASH.mat corresponds to the desired set.
    spec_ids_hash = string2hash(strjoin(spec_ids, '_'));
    
    fps_mat_fn_basename = sprintf('fps_%0.f.mat', spec_ids_hash);
    fps_mat_fn = fullfile(mat_fps_dir, fps_mat_fn_basename);
    if exist(fps_mat_fn, 'file')  
        if verbose 
            fprintf('Load "%s" file.\n', fps_mat_fn_basename);
        end % if
        
        load(fps_mat_fn, 'fps_struct');
    else
        n_spec = length(spec_ids);
        
        if verbose
            fprintf('Load separated *.fpt files: # = %d. Create "%s" subsequently.\n', n_spec, ...
                fps_mat_fn_basename);
        end % if 
        
        % Load the fingerprint mask and determine fingerprint dimension
        mask = load_fingerprint_mask(fullfile(raw_fps_dir, "fingerprints.mask"));
        n_fps = length(mask);
        
        fps_struct = struct('fps', false(n_fps, n_spec), 'mask', mask);
        
        for ii = 1:n_spec            
            fps_spec_fn = fullfile(raw_fps_dir, strcat(spec_ids{ii}, '.fpt'));
            fid = fopen(fps_spec_fn, 'r');
            if fid == -1
                error('loadFingerprints:NoSuchFileOrDirectory', 'Cannot open fps-file "%s".', ...
                    fps_spec_fn);
            end % if
            
            % parse fingerprint vector 
            fps_struct.fps((fgetl(fid) - '0') == 1, ii) = true;
            
            fclose(fid);
            
            if verbose && ((mod(ii, 500) == 0) || (ii == n_spec))
                fprintf('Processed spectra: %d/%d.\n', ii, n_spec);
            end % if
        end % if
        
        fps = fps_struct.fps;
        mask = fps_struct.mask;
        
        % save the fingerprint structure fur future use
        save(fps_mat_fn, 'fps_struct', '-v7.3');
    end % if
    
    if nargout == 1
        fps = fps_struct;
    elseif nargout == 2
        fps = fps_struct.fps;
        mask = fps_struct.mask;
    end % if
end % function
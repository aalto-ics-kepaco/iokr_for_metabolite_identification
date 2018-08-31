function cmps = loadCompoundList(input_dir, format)
%% LOADCOMPOUNDLIST Loads the csv-file containing the list of ids and cmps.
    switch format
        case 'csi_fingerid'
            fn = input_dir + '/compounds';
            if ~ exist(fn, 'file') 
                error('loadCompoundList:FileNotFound', ...
                    ['Cannot find training ids file: %s. Please run "fingerID list-compounds > compounds"' ...
                    'or "fingerID list-independent-compounds --indep-set=SET > SET/compounds" first.'], fn);
            end % if
            
            fid = fopen(fn, 'r');
            if (~ fid)
                error('loadCompoundList:IOError', 'Cannot open "%s",', fn);
            end % if

            cmps = textscan(fid, '%s%s%s', 'Delimiter', '\t');
            cmps = table(cmps{1}, cmps{2}, cmps{3}, 'VariableNames', {'spec_id', 'inchikey2D', 'inchi2D'});
                        
            fclose(fid);
        otherwise
            error('loadCompoundList:InvalidArgument', '"%s" is not a valid format.', format)
    end % swtich
end % functions,
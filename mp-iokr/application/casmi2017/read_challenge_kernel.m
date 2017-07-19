function [KX, KX_names] = read_challenge_kernel (filename)
    tic;

    fid = fopen (filename, 'r');
    if (fid < 0)
        error ('Cannot open file: %s', filename);
    end % if
    
    KX_names = cell(0);
    KX       = cell(0);
    k        = 1;
    
    line = fgetl (fid);
    while (ischar (line))
        line_splitted = strsplit (line, '\t');
        
        KX_names{k} = line_splitted{1};
        KX{k}       = cellfun (@(x) str2double (x), line_splitted(2:end))';
        
        line = fgetl(fid);
        k    = k + 1;
    end % if
    
    fclose (fid);
    
    toc;
end % function 

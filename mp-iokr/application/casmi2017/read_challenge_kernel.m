function [KX, KX_names] = read_challenge_kernel(filename, split_char, verbose)
    %% READ_CHALLENGE_KERNEL reads in kernels line by line
    % Function to read training-vs-test or test vs test kernels calculated
    % with the CSI:FingerID framework.
    %
    % Structure of the training-vs-test kernel files:
    %   KERNEL_1[split_char]k_1(x_1, x')[split_char]k_1(x_2, x')...[split_char]k_1(x_n, x')\n
    %   KERNEL_2[split_char]k_2(x_1, x')[split_char]k_2(x_2, x')...[split_char]k_2(x_n, x')\n
    %   ...
    %
    %   with n being the number of training examples and x' being the
    %   specific test example.
    %
    %   E.g.
    %   PPKr 0.12 0.23 ... 0.31
    %   MLIP 0.52 0.4233 ... 0.769
    %   ...
    if nargin < 2
        split_char = '\t';
    end % if
    
    if nargin < 3
        verbose = true;
    end % if

    if verbose
        tic;
    end % if

    fid = fopen(filename, 'r');
    if (fid < 0)
        error('read_challenge_kernel:NoSuchFileOrDirectory', 'Cannot open file: %s', filename);
    end % if
    
    KX_names = cell(0);
    KX       = cell(0);
    k        = 1;
    
    line = fgetl(fid);
    while (ischar(line))
        line_splitted = strsplit(line, split_char);
        
        KX_names{k} = line_splitted{1};
        KX{k}       = cellfun (@(x) str2double(x), line_splitted(2:end))';
        
        line = fgetl(fid);
        k    = k + 1;
    end % if
    
    fclose(fid);
    
    if verbose
        toc;
    end % if
end % function 

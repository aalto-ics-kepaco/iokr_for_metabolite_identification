%% SAVEKERNEL Function to save a Kernel matrix into a file
% This function saves a given kernel matrix and it's header (if wanted)
% into a *.txt or *.mat file. The output-file format depends on the
% file-ending your provie with the filename.
%
% Input:
%   oFilename   ... string containing the path to the kernel-matrix file
%   KX          ... kernel-matrix
%   header[={}] ... cell array containing the name of each example in the
%                   kernel matrix
%                   NOTE: If the output file is *.txt, then the first row
%                   the text file will be the header (starting with '#'). 
%                   If the output file is *.mat, then header will be stored
%                   as 'header' in the output-struct.
function saveKernel (oFilename, KX, header)
    tic;
    
    if (nargin < 3)
        header = [];
    end % if
    
    [~, ~, fileExt] = fileparts (oFilename);
    switch (fileExt)
        case '.txt'    
            if (~isempty (header))
                
                fid = fopen (oFilename, 'w');
                if (~fid)
                    error ('saveKernel:InvalidInput', 'Can not open output-file: %s.\n', oFilename);
                end % if

                fprintf (fid, '# ');
                fprintf (fid, '%s ', header{1:(end-1)});
                fprintf (fid, '%s', header{end});
                fprintf (fid, '\n');

                fclose (fid);

                dlmwrite (oFilename, KX, 'delimiter', ' ', '-append', 'precision', '%.18f');    
            else 
                save (oFilename, 'KX', '-ascii', '-double');    
            end % if
        case '.mat'
            save (oFilename, 'KX', '-v7.3');
            if (~isempty (header))
                save (oFilename, 'header', '-append');
            end % if
        otherwise 
            error ('saveKernel:InvalidInput', ...
                'The kernel-matrix files is not valid. Only .txt and .mat are allowed.\n');
    end % switch
   
    fprintf ('Saving the Kernel-matrix took %.3fs.\n', toc);
end % function
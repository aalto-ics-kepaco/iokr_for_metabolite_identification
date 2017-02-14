function saveKernel (filename, KX, header, overwrite)
%% SAVEKERNEL Function to save a Kernel matrix into a file
%   SAVEKERNEL (FILENAME, KX) saves the kernel-matrix KX into the file at
%   FILENAME. The provided file-extention of FILENAME determines the
%   file-format used to store the kernel-matrix (see 'notes'). 
%
%   SAVEKERNEL (..., HEADER) Additionally a cell-array HEADER can be
%   provided, e.g. containing the identifiers of the examples in the
%   kernel-matrix. Depending on the output file-format the HEADER will be
%   handled in a different way (see 'notes').
%
%   SAVEKERNEL (..., HEADER, OVERWRITE) If the binary value OVERWRITE is
%   set to true than a potentially already existing kernel-matrix file will
%   be overwritten _without_ any additional message.
%
%   INPUTS:
%       filename        Filename to store the kernel-matrix.
%       KX              M x N kernel-matrix
%       header                                                         [{}]
%                       Cell-array containing the elements which will be
%                       put into the header line. If 'isempty (header)' is
%                       true. _No_ header will be written out.
%       overwrite                                                   [false]
%                       Binary value indicating whether a possibly already
%                       existing kernel-matrix file will be overwritten.
%
%   NOTES:
%       Depending on the file extention (either '.txt' or '.mat') the
%       kernel-matrix will be stored in a different way:
%           * .txt:     If a header is provided it will be written in the
%                       first line of the kernel-matrix file. The header
%                       line begins with a '#' _followed_ by a space. The
%                       elements of the header are written out separated by 
%                       'space'. 
%                       The kernel-matrix is written out line-by-line using
%                       'space' as separator between the different values.
%                       The precision of the values is '%1.16e'. 
%           * .mat:     The kernel-matrix and if provided the header are
%                       stored using matlabs SAVE function. The
%                       kernel-matrix has the identifier 'KX' and the
%                       header 'header'.
%
%   See also LOADKERNEL.
    tic;
    
    %% Check the input and set arguments to defaults if needed
    if (nargin < 3) ; header    = {}    ; end % if
    if (nargin < 4) ; overwrite = false ; end % if
        
    %% Store the provided kernel-matrix
    if (exist (filename, 'file') && ~ overwrite)
        warning ('%s does already exist. Nothing is written out.', filename);
    else
        [~, ~, fileExt] = fileparts (filename);
        switch (fileExt)
            case '.txt'    
                if (~ isempty (header))
                    % Write out header
                    fid = fopen (filename, 'w');
                    if (~ fid)
                        error ('saveKernel:InvalidInput', 'Can not open output-file: %s.\n', filename);
                    end % if

                    fprintf (fid, '# ');
                    fprintf (fid, '%s ', num2str (header{1:(end-1)}));
                    fprintf (fid, '%s', num2str (header{end}));
                    fprintf (fid, '\n');

                    fclose (fid);
                    
                    % Write out data
                    dlmwrite (filename, KX, '-append', 'delimiter', ' ', 'precision', '%1.16e');    
                else
                    % Write out data
                    dlmwrite (filename, KX, 'delimiter', ' ', 'precision', '%1.16e');    
                end % if
            case '.mat'
                save (filename, 'KX', '-v7.3');
                if (~ isempty (header))
                    save (filename, 'header', '-append', '-v7.3');
                end % if
            otherwise 
                error ('saveKernel:InvalidInput', ...
                    'The kernel-matrix file is not valid. Only *.txt and *.mat are allowed as file-extensions.\n');
        end % switch
    end % if 

    fprintf ('Saving the Kernel-matrix took %.3fs.\n', toc);
end % function
%% LOADKERNEL Function to load a kernel-matrix
% This function can handle different kernel-matrix files depending on the 
% kernel-file's extension: 
% 1) ... '.txt' (e.g. output of Kai's tool)
%    The file is loaded using 'readtable' (probably line-by-line and very
%    slow). The file can contain a header line ('# mol1 mol2 ...'). 
% 2) ... '.mat' (some output from matlab)
%    The file is loaded using 'load' and should contain at least the kernel
%    matrix 'KX'. It can contain an element 'header' too.
%
% Input:
%   filename       ... string containing the path to the kernel-matrix file
%   hasHeader[=1]  ... binary indicating whether the kernel-file contains a
%                      header or not
%   loadOnlyHeader ... binary indicating whether only the header of the
%                      kernel-file should be loaded
%   [=0]
% Output:
%   kernel   ... kernel-matrix
%   [header] ... header line of the kernel matrix as cell-array
function [KX, header] = loadKernel (filepath, hasHeader, loadOnlyHeader)
    tic;
    if (~ exist (filepath, 'file'))
        error ('loadKernel:FileNotFound', ...
            'Kernel-file: %s does not exist\n.', filepath);
    end % if
    
    if (nargin < 2) ; hasHeader = 1 ; end % if
    if (nargin < 3) ; loadOnlyHeader = 0; end % if
    if ((~ hasHeader) && (nargout > 1))
        warning ('You indicate that the kernel-file does not contain a header, but you request it. So an empty header will be returned.');
        header = {};
    end % if
    
    [~, ~, fileExt] = fileparts (filepath);
    switch (fileExt)
        case '.txt'
            % Read kernel-matrix as table (e.g: from Kai's tool)
            if (loadOnlyHeader)
                KX = [];
            else
                KX = readtable (filepath, ...
                    'ReadVariableNames', 0, 'HeaderLines', hasHeader, 'Delimiter', 'space');
                KX = table2array (KX);
            end % if

            if (hasHeader && (nargout > 1))
                % If header should be reported, read the first line of the
                % kernel-matrix file and parse it
                fid = fopen (filepath, 'r');

                header = fgetl (fid);
                header = strsplit (header);
                header = header(2:end);

                fclose (fid);
            end % if
        case '.mat'
            % Load header if needed
            if (hasHeader && (nargout > 1))
                load (filepath, 'header');
                if (~ exist ('header', 'var')) 
                    error ('loadKernel:InvalidInput', ...
                        'No header in the .mat-file. Please store the header using "header" as variable name.\n');
                end % if
            end % if
            
            % Read kernel-matrix as matlab-object
            if (loadOnlyHeader)
                KX = [];
            else
                load (filepath, 'KX');
                if (~ exist ('KX', 'var'))
                    error ('loadKernel:InvalidInput', ...
                        'No kernel-matrix in the .mat-file. Please store the kernel-matrix using "KX" as variable name.\n');
                end % if
            end % if
        otherwise
            error ('loadKernel:InvalidInput', ...
                'The kernel-matrix files is not valid. Only .txt and .mat are allowed.\n');
    end % switch
      
    fprintf ('Loading the kernel took %.3fs.\n', toc);
end % function
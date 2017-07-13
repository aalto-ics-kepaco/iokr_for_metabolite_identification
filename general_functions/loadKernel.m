function [KX, header] = loadKernel (filepath, loadHeaderOnly)
%% LOADKERNEL Function to load a kernel-matrix
%   [KERNELMATRIX] = loadKernel (FILEPATH) loads the kernel matrix stored 
%   at FILEPATH. 
%
%   [KERNELMATRIX] = loadKernel (FILEPATH, LOADHEADERONLY) by setting the 
%   binary value LOADHEADERONLY to true _only_ the header of the 
%   kernel-matrix is loaded from the file. An empty KERNELMATRIX is returned.
%
%   [KERNELMATRIX, HEADER] = loadKernel (...) If the kernel-matrix file 
%   contains a headerm e.g. some kind of identifier for each example in the 
%   kernel-matrix. It is returned and stored into HEADER. See 'notes' for 
%   more information about the encoding of the header. 
%
%   INPUTS:
%       filepath        Path to the kernel-matrix file. See also 'notes'.
%       loadHeaderOnly                                              [false]
%                       Binary indicating whether only the header of the
%                       kernel-matrix, e.g. containing the examples'
%                       identifier, should be loaded. 
%
%   OUTPUTS:
%       kernelMatrix    M x N matrix
%       header          Cell-array containing the examples' identifiers.
%
%   NOTES:
%       This function can handle different kernel-matrix files depending on the 
%       kernel-file's extension: 
%       1) ... '.txt' (e.g. output of Kai's tool)
%              The file is loaded using 'readtable' (probably line-by-line 
%              and very slow). The file can contain a header line, e.g.
%              '# mol1 mol2 ...', which start with a '#' whereby the different
%              identifiers are separated by SPACE. 
%       2) ... '.mat' (some output from matlab)
%              The file is loaded using 'load' and should contain at least 
%              the kernel matrix 'KX'. It can contain an element 'header', 
%              which is returned if the HEADER is requested.
%
%       The function SAVEKERNEL can be used to store kernel in the propper
%       format. 
%
%   See also SAVEKERNEL.
    tic;
    if (~ exist (filepath, 'file'))
        error ('loadKernel:InvalidInput', ...
            'Kernel-matrix file: %s does not exist.\n', filepath);
    end % if
    
    if (nargin < 2) ; loadHeaderOnly = false; end % if
    header = {};
    
    [~, ~, fileExt] = fileparts (filepath);
    switch (fileExt)
        case '.txt'
            % Check whether the kernel-matrix file contains a header. For
            % that read the first line of the kernel-matrix file and parse
            % it. If the line starts with a '#' it will be interpreted as
            % header line.
            fid = fopen (filepath, 'r');
            if (~ fid)
                error ('loadKernel:InvalidInput', ...
                    'Cannot open kernel-matrix file: %s.\n', filepath);
            end %if

            firstLine    = fgetl (fid);
            firstLine    = strsplit (firstLine);
            hasHeader    = strcmp (firstLine(1), '#');
            
            % Load the header if requested
            if (nargout > 1)          
                if (hasHeader)
                    header = firstLine(2:end); 
                else
                    warning ('Header requested but not provided. An empty header will be returned.');
                end % if
            end % if        
            
            fclose (fid);
                
            % Read kernel-matrix as table, e.g the output of Kai's tool
            if (loadHeaderOnly)
                KX = [];
            else
                KX = readtable (filepath, ...
                    'ReadVariableNames', 0, 'HeaderLines', double (hasHeader), 'Delimiter', 'space');
                KX = table2array (KX);
            end % if
        case '.mat'
            % Load header if requested
            if (nargout > 1)
                load (filepath, 'header');
                
                % Check whether the kernel-matrix file contains a header
                if (~ exist ('header', 'var'))
                    warning ('Header requested but not provided. An empty header will be returned.');
                end % if
            end % if
            
            % Read kernel-matrix as matlab-object
            if (loadHeaderOnly)
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
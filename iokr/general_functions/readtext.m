function [ t ] = readtext( filename )
%======================================================
% DESCRIPTION:
% Read data from a file containing one string on each line
%
% INPUTS:
% filename:     string corresponding to the file name
%
% OUTPUTS:
% t:            cell of strings
%
%======================================================

    fileID = fopen(filename);
    C = textscan(fileID, '%s'); 
    t = C{1};
    fclose(fileID);
    
end

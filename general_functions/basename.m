function filename = basename (fullFn)
%% BASENAME removes the filepath and -extention of a given filename.
%   filename = BASENAME (fullFn) returns the filename extracted from
%   fullFn. 
%   Example:
%       basename ('/hello/world.txt') returns 'world'
    [~, filename, ~] = fileparts (fullFn); 
end % function
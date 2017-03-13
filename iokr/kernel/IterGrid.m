function [ params ] = IterGrid( ky_param )
%======================================================
% DESCRIPTION:
% Generates all the possible parameter combinations for the kernel parameters
%
% INPUTS:
% ky_param:     struct array of size 1*1 containing the kernel options and the 
%               parameter values to consider for the parameter selection 
%
% OUTPUTS:
% params:       struct array of size 1*n where n is the number of possible
%               parameter combinations. Each value corresponds to a parameter combination.
%
% EXAMPLE:
%
%   >> ky_param = struct('type','gaussian','base_kernel','linear', 'gamma', [1 2 4 8]);
%   >> params = IterGrid( ky_param );
%
% params is a 4*1 struct array with fields: type, base_kernel and gamma.
%
% params(1) corresponds to struct('type','gaussian','base_kernel','linear', 'gamma', 1)
% params(2) corresponds to struct('type','gaussian','base_kernel','linear', 'gamma', 2)
% params(3) corresponds to struct('type','gaussian','base_kernel','linear', 'gamma', 4)
% params(4) corresponds to struct('type','gaussian','base_kernel','linear', 'gamma', 8)
%
%======================================================
    

    fields = fieldnames(ky_param); % fields of the struct array
    num_fields = length(fields);
     
    % extract the values for each field in a cell array
    ky_values = cell(num_fields,1);
    for i = 1:num_fields
        
        field_value = ky_param.(fields{i});
        
        if ischar(field_value)
            ky_values{i} = {field_value};
        else
            ky_values{i} = num2cell(field_value);
        end
        
    end

    % compute all the possible pairs of parameters
    ix = cellfun(@(x) 1:numel(x), ky_values, 'UniformOutput', false);

    [ix{:}] = ndgrid(ix{:});

    A = cell(numel(ix{1}), num_fields);
    for i = 1:num_fields

        A(:,i) = reshape(ky_values{i}(ix{i}),[],1);

    end

    % convert the cell array to a struct array
    params = cell2struct(A, fields, 2);
    
end

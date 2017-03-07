function [ val_lambda, params ] = IterGrid( param_grid, KY_opt )
%======================================================
% DESCRIPTION:
% Generates all the possible parameter combinations for the parameter selection
%
% INPUTS:
% param_grid:   struct array of size 1*1 indicating the parameter grid to explore 
%
% OUTPUTS:
% params:       struct array of size 1*n where n is the number of possible
%               parameter combinations. Each value corresponds to a parameter combination.
%
% EXAMPLE:
%   >> param_grid = struct('lambda', [1 10], 'gamma', [2 3 4]);
%   >> t = IterGrid(param_grid)
%   ans = 
% 
%   1x6 struct array with fields:
% 
%     lambda
%     gamma
%
%   >> t(1)
%   ans = 
% 
%     lambda: 1
%      gamma: 2
%
%   >> t(2)
% 
%   ans = 
% 
%     lambda: 10
%      gamma: 2
%
%======================================================

    % remove the regularization parameter
    val_lambda = param_grid.lambda;
    param_grid = rmfield(param_grid,'lambda');

    f = fieldnames(param_grid);

    if isempty(f)
        % in the case of kernel without parameters, params is a 1*1 struct
        % array containing the kernel type and eventually the base kernel
        
        params = KY_opt;
        
    else
        
        f_ky = fieldnames(KY_opt);
        f2 = {f_ky{:},f{:}}';
        % otherwise we enumerate all the possible combinations for the
        % kernel parameter(s)
        
        % Extracts the values for the different fields in a cell array
        c = cellfun(@(x) getfield(param_grid,x), f,'UniformOutput', false);
        c_ky = cellfun(@(x) getfield(KY_opt,x), f_ky,'UniformOutput', false);

        % Generates the parameter grid
        b = cell(numel(c),1);
        [b{:}] = ndgrid(c{:});
        n_param = length(b{1}(:));

        % Builds a struct array containing all the possible parameter combinations
        
        params = cell2struct(cell(length(f2),1), f2, 1);
        
        for k = 1:n_param
            v = cellfun(@(x) x(k), b,'UniformOutput', false);
            v2 = {c_ky{:},v{:}}';
            params(k) = cell2struct(v2, f2, 1);
            
        end
    end
    
end
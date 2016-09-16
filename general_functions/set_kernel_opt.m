function [ kernel_opt ] = set_kernel_opt( kernel_type, kernel_param )
%======================================================
% DESCRIPTION:
% Set the kernel_opt structure containing the kernel type and parameter(s)
%
% INPUTS:
% kernel_type:     kernel type ('gaussian', 'polynomial' or 'linear'
% kernel_param:    kernel parameter(s)
%       linear kernel: no parameter needed
%       polynomial kernel: vector of two values, the first one corresponds 
%                          to the offset and the second to the degree
%       gaussian kernel: gamma parameter
%
% OUTPUTS:
% kernel_opt: structure containing the kernel type and parameter(s)
%
%======================================================

    kernel_opt.type = kernel_type;
    switch kernel_opt.type
        case 'gaussian'
            kernel_opt.gamma = kernel_param;
            
        case 'polynomial'
            kernel_opt.scale = 1;
            kernel_opt.offset = kernel_param(1);
            kernel_opt.degree = kernel_param(2);
    end

end
    
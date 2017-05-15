function [ Y_norm ] = norma( Y, mean_Y, center_opt )
%======================================================
% DESCRIPTION:
% Vector centering and normalization
%
% INPUTS:
% Y:                matrix of size d*n, in which each column corresponds to 
%                   a vector yi (Y = [y1,...,yn])
% mean_Y:           vector containing the mean used for centering
% center_opt:       binary value indicating if the vectors should be
%                   centered (1) or not (0)
%
% OUTPUTS:
% Y_norm:           matrix containing the centered and normalized vectors
%
%====================================================== 

    [d,n] = size(Y);
    
    if center_opt ==1
        Y_c = Y - repmat(mean_Y,1,n);
    else
        Y_c = Y;
    end
    
    D = sqrt(sum(Y_c.*Y_c));
    Y_norm = Y_c ./ repmat(D,d,1); 

end


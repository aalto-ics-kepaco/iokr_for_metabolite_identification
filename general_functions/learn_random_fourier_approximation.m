function W = learn_random_fourier_approximation (Y, gamma, D)
% Random Fourier features: approximates the feature vectors of a Gaussian
% kernel with parameter gamma.
% The dimension of the approximated feature vectors will be 2*D.
% Higher is D, better will be the approximation.
%
% Reference paper: Rahimi, A. & Recht, B. Random features for large-scale
% kernel machines. NIPS, 2007.

    d = size(Y,1); % number of molecular properties
    
    % Draw a random matrix from the standard normal distribution
    W = sqrt(2 * gamma) * randn(D, d);
end
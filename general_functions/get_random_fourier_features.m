function Psi = get_random_fourier_features (Y, W)
    D = size (W, 1);

    % Approximated feature vectors
    Psi = sqrt(1/D) * [cos(W*Y); sin(W*Y)];
end % function 
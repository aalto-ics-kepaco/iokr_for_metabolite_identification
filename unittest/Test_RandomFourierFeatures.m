classdef Test_RandomFourierFeatures < matlab.unittest.TestCase
    methods (Test)
        function testApproximationOnSmallDataset (tc)
            rng (10);
            
            close all;
            
            % Load fingerprints
            inputDir = '/home/bach/Documents/studies/doctoral/data/csi_fingerid/example/iokr/input/';
            load ([inputDir '/compound_info.mat'], 'dt_inchi_mf_fp');
            
            Y      = full (dt_inchi_mf_fp.fp_masked)';
            [d, n] = size (Y);
            clear dt_inchi_mf_fp;
            
            % Estimate the optimal gamma
            ky_param       = struct ('type', 'gaussian', 'base_kernel', 'linear');
            ky_param.gamma = select_gamma_entropy (Y, ky_param);
            
            % Calculate the true kernel
            K_orig = build_kernel (Y, Y, ky_param);
            
            % Get the random fourier features for different D values
            val_D  = [10, 100, 500, 1e3, 5e3, 1e4, 5e4, 1e5] * 2;
            err    = zeros (1, length (val_D));
            D_appr = zeros (1, length (val_D));
            
            figure;
            subplot (2, 5, 1);
            imagesc (K_orig); title ('Original');
            
            for ii = 1:length (val_D)
                rff = RandomFourierFeatures (d, val_D(ii), ky_param.gamma);
                Psi = rff.getRandomFourierFeatures (Y);
                
                K_appr     = build_kernel (Psi, Psi, struct ('type', 'linear'));
                err(ii)    = sum (sum (abs (K_appr - K_orig))) / (n * n);
                D_appr(ii) = size (Psi, 1);
                subplot (2, 5, ii + 1);
                imagesc (K_appr); title (sprintf ('Approx: D = %d', D_appr(ii)));
            end % for
            
            figure;
            semilogx (D_appr, err, '--*'); 
            xlabel ('D'); ylabel ('Error'); grid;
        end % function 
    end % Test methods
end % class

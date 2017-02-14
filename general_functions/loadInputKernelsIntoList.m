function [KX_list, param] = loadInputKernelsIntoList (inputDir, param)
%% LOADINPUTKERNELLIST
%   Function to load a list of kernels. Either a single input kernel or a 

    switch (upper (param.data_param.inputKernel))
        case upper (param.data_param.availInputKernels)
            disp (['Evaluation using a single kernel: ', param.data_param.inputKernel]);

            KX_list = { loadKernel(strcat (inputDir, '/', upper (param.data_param.inputKernel), '.mat')) };

            % If a single kernel is used, than we force the "mkl" option to
            % be 'unimkl'. The weight for the kernel will be 1.
            param.mp_iokr_param.mkl = 'unimkl';
        case upper ({'unimkl', 'alignf'})
            disp (['Evaluation using multiple kernel learning: ', param.data_param.inputKernel]);

            KX_list = cellfun (@(kernelName) loadKernel (strcat (inputDir, '/', kernelName, '.mat')), ...
                upper (param.data_param.availInputKernels), 'UniformOutput', false);

            param.mp_iokr_param.mkl = lower (param.data_param.inputKernel);
        otherwise
            error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', ...
                '%s: Not a valid input kernel. See MP_IOKR_DEFAULTS for a list of available kernels.', ...
                param.data_param.inputKernel);
    end % switch
end % function

%             warning ('Configured to read .txt kernel-files!');
%             KX_list = { loadKernel(strcat (inputDir, '/input_kernels/', upper (param.data_param.inputKernel), '.txt'), 0) };

%             KX_list = cellfun (@(kernelName) loadKernel (strcat (inputDir, '/input_kernels/', kernelName, '.txt'), 0), ...
%                 upper (param.data_param.availInputKernels), 'UniformOutput', false);
%             warning ('Configured to read .txt kernel-files!');
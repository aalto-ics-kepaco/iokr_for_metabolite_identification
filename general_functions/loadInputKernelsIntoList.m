function [KX_list, param] = loadInputKernelsIntoList (inputDir, param, fileExtension)
%% LOADINPUTKERNELSINTOLIST Loads a list of kernels matrices.
    if nargin < 3
        fileExtension = '.mat';
    end % if

    switch (upper (param.data_param.inputKernel))
        case upper (param.data_param.availInputKernels)
            disp (['Evaluation using a single kernel: ', param.data_param.inputKernel]);

            KX_list = { loadKernel(strcat (inputDir, '/', upper (param.data_param.inputKernel), fileExtension)) };

            % If a single kernel is used, than we force the "mkl" option to
            % be 'unimkl'. The weight for the kernel will be 1.
            param.mp_iokr_param.mkl = 'unimkl';
        case upper ({'unimkl', 'alignf'})
            disp (['Evaluation using multiple kernel learning: ', param.data_param.inputKernel]);

            KX_list = cellfun (@(kernelName) loadKernel (strcat (inputDir, '/', kernelName, fileExtension)), ...
                upper (param.data_param.availInputKernels), 'UniformOutput', false);

            param.mp_iokr_param.mkl = lower (param.data_param.inputKernel);
        otherwise
            error ('IOKR_MP_reverse_feat_evaluation:InvalidInput', ...
                '%s: Not a valid input kernel. See MP_IOKR_DEFAULTS for a list of available kernels.', ...
                param.data_param.inputKernel);
    end % switch
end % function
function convert_kernel_files_to_mat (wdir)
    %% CONVERT_KERNEL_FILES_TO_MAT Load kernels as 'txt' and save it as 'mat'
    % All *.txt files in 'kernel_dir' are interpreted as kernel matrices
    % calculated using csi_fingerid. The matrices are loaded and saved as
    % mat-files.
    
    kernel_files = dir ([wdir, '/*.txt']);
    for kernel_file = kernel_files'
        fprintf ('Load kernel: %s\n', basename (kernel_file.name));
        [KX, header] = loadKernel ([wdir, '/', kernel_file.name]);
        saveKernel ([wdir, '/', basename(kernel_file.name), '.mat'], ...
            KX, header);
    end % for
end % function
function convert_kernel_files_to_mat (kernel_dir)
    %% CONVERT_KERNEL_FILES_TO_MAT Load kernels as 'txt' and save it as 'mat'
    % All *.txt files in 'kernel_dir' are interpreted as kernel matrices
    % calculated using csi_fingerid. The matrices are loaded and saved as
    % mat-files.
    kernel_files = dir ([kernel_dir, '/*.txt']);
    for kernel_file = kernel_files'
        [KX, header] = loadKernel ([kernel_dir, '/', kernel_file.name]);
        saveKernel ([kernel_dir, '/', basename(kernel_file.name), '.mat'], ...
            KX, header);
    end % for
end % function
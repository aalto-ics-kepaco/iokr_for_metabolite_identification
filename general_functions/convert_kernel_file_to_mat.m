function convert_kernel_file_to_mat (wdir, kernel_name)
    %% CONVERT_KERNEL_FILES_TO_MAT Load kernels as 'txt' and save it as 'mat'
    % All *.txt files in 'kernel_dir' are interpreted as kernel matrices
    % calculated using csi_fingerid. The matrices are loaded and saved as
    % mat-files.
    
    fprintf ('Load kernel: %s\n', kernel_name);
    [KX, header] = loadKernel ([wdir, '/', kernel_name, '.txt']);
    saveKernel ([wdir, '/', kernel_name, '.mat'], KX, header);
end % function
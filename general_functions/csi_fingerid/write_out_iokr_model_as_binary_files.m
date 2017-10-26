function write_out_iokr_model_as_binary_files (train_model, kernel_names, odir)
    assert (length (train_model.process_input) == length (kernel_names))

    %% Process (input) row-means ,diag_c and mkl-weights
    for i_kernel__ = 1:length (kernel_names)
        fid__ = fopen ( ...
            [odir, '/KX_mean_', basename(kernel_names{i_kernel__}), '.bin'], 'w');
        assert (fid__ > 0);
        fwrite (fid__, train_model.process_input(i_kernel__).mean, 'double', 'ieee-le');
        fclose (fid__);

        fid__ = fopen ( ...
            [odir, '/KX_diag_c_', basename(kernel_names{i_kernel__}), '.bin'], 'w');
        assert (fid__ > 0);
        fwrite (fid__, train_model.process_input(i_kernel__).diag_c, 'double', 'ieee-le');
        fclose (fid__);

        fid__ = fopen ( ...
            [odir, '/mkl_w_', basename(kernel_names{i_kernel__}), '.bin'], 'w');
        assert (fid__ > 0);
        fwrite (fid__, train_model.process_input(i_kernel__).w, 'double', 'ieee-le');
        fclose (fid__);    
    end % for

    %% Process triangular matrix
    fid__ = fopen ([odir, '/C.bin'], 'w');
    assert (fid__ > 0);
    Ct__ = train_model.C';
    fwrite (fid__, Ct__(Ct__ ~=0), 'double', 'ieee-le');
    fclose (fid__);

    %% Process row-means (input)
    fid__ = fopen ([odir, '/KY_mean.bin'], 'w');
    assert (fid__ > 0)
    fwrite (fid__, train_model.process_output.mean, 'double', 'ieee-le');
    fclose (fid__);

    fid__ = fopen ([odir, '/KY_diag_c.bin'], 'w');
    assert (fid__ > 0);
    fwrite (fid__, train_model.process_output.diag_c, 'double', 'ieee-le');
    fclose (fid__);

    %% Gamma
    fid__ = fopen ([odir, '/KY_gamma.bin'], 'w');
    assert (fid__ > 0)
    fwrite (fid__, train_model.KY_par.gamma, 'double', 'ieee-le');
    fclose (fid__);
end % function
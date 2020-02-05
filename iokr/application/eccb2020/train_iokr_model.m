function [ ] = train_iokr_model (inputDir, outputDir, ionization_mode)
%======================================================
% DESCRIPTION:
% Script for running IOKR on the
%
% INPUTS:
% inputDir:        string, base-dictionary containing the input data
% result_dir:      string, base-dictionary where the output should be stored
% ionization_mode: string, which ionization to predict the scores for {'positive', 'negative'}
%
%======================================================

    %--------------------------------------------------------------
    % Set up parameters
    %--------------------------------------------------------------
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
        {'debug_param', 'opt_param', 'iokr_param', 'data_param', 'ky_param'});
    
    % Number of cross-validation folds for the hyper-parameter
    % optimization, e.g. finding lambda (regularization parameter).
    param.opt_param.nInnerFolds = 10;

    % Output kernel: kappa_y
    param.ky_param.representation  = 'kernel';
    param.ky_param.type            = 'gaussian';
    param.ky_param.base_kernel     = 'tanimoto';
    param.ky_param.param_selection = 'entropy';
    
    param.debug_param.verbose = true;
    
    %--------------------------------------------------------------
    % Load and prepare data
    %--------------------------------------------------------------
    % Information about MS/MS spectra ids, CCS values and molecular
    % formulas for the training compounds. The order of the comounds in
    % this table, corresponds to the order of the rows and columns in the
    % input kernel matrices. 
    cmps = readtable(fullfile(inputDir, 'compounds.csv'));
    
    % Fingerprints
    fps_dir = fullfile(inputDir, 'fingerprints');
    [Y_train, ~] = loadFingerprints(fps_dir, fps_dir, cmps.spec_id);  % Y with shape = (n_fps, n_samples)
         
    % Input kernels
    kernel_dir = fullfile(inputDir, 'kernels');
    kernel_files = dir(fullfile(kernel_dir, '*.mat'));
    param.data_param.availInputKernels = arrayfun(@(file) basename(file.name), ...
        kernel_files, 'UniformOutput', false);
    KX_list_train = loadInputKernelsIntoList(kernel_dir ,param, '.mat');
    
    %--------------------------------------------------------------
    % Train IOKR Model
    %--------------------------------------------------------------
    iokr_model = Train_IOKR (KX_list_train, Y_train, ...
        param.ky_param, param.opt_param, param.iokr_param, param.debug_param.verbose);
    
    save ()
end
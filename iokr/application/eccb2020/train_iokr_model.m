function [ ] = train_iokr_model (inputDir, outputDir)
%======================================================
% DESCRIPTION:
% Script to train an IOKR model
%
% INPUTS:
% inputDir:       string, base-dictionary containing the input data
% outputDir:      string, base-dictionary where the output should be stored
%                 (Default is the same as the 'inputDir')
%
% Folders and files that need to exists for this script to run:
%   inputDir/kernels:       Directory containing the training kernels (*.mat)
%   inputDir/compounds.csv: Output if 'fingerID list-compounds', a
%                           tab-separated header-less list with three columns
%                           <SpecID>\t<InChIKey1>\t<InChI2D> 
%   inputDir/fingerprints:  Directory containing the fingerprint files (*.fps)
%                           for all training compounds. There must be one file
%                           per <SpecID>.
%======================================================
    if nargin < 2
        outputDir = inputDir;
    end % if

    %--------------------------------------------------------------
    % Set up parameters
    %--------------------------------------------------------------
    param = MP_IOKR_Defaults.setDefaultsIfNeeded (struct(), ...
        {'debug_param', 'opt_param', 'iokr_param', 'data_param', 'ky_param'});
    
    % Number of cross-validation folds for the hyper-parameter
    % optimization, e.g. finding lambda (regularization parameter).

    % Output kernel: kappa_y
    param.ky_param.representation  = 'kernel';
    param.ky_param.type            = 'gaussian';
    param.ky_param.base_kernel     = 'tanimoto';
    param.ky_param.param_selection = 'entropy';
    
    param.iokr_param.model_representation = 'Chol_decomp_of_C';
    
    param.debug_param.verbose = true;
    
    %--------------------------------------------------------------
    % Load and prepare data
    %--------------------------------------------------------------
    cmps = readtable(fullfile(inputDir, 'compounds.csv'), 'Delimiter', '\t', 'ReadVariableNames', false);
    
    % Fingerprints
    fps_dir = fullfile(inputDir, 'fingerprints');
    [Y_train, ~] = loadFingerprints(fps_dir, fps_dir, cmps.Var1);  % Y with shape = (n_fps, n_samples)
    [n_fps, n_samples] = size(Y_train);
    
    fprintf('Number of fingerprint vectors: %d\n', n_samples);
    fprintf('Fingerprint dimension: %d\n', n_fps);
         
    % Input kernels
    kernel_dir = fullfile (inputDir, 'kernels');
    kernel_files = dir (fullfile (kernel_dir, '*.mat'));
    param.data_param.availInputKernels = arrayfun(@(file) basename(file.name), ...
        kernel_files, 'UniformOutput', false);
    KX_list_train = loadInputKernelsIntoList(kernel_dir ,param, '.mat');
    
    [n_samples_kernel, ~] = size(KX_list_train{1});
    fprintf('Number of rows in training kernel matrix: %d\n', n_samples_kernel);
    
    %--------------------------------------------------------------
    % Train IOKR Model
    %--------------------------------------------------------------
    model = Train_IOKR (KX_list_train, Y_train, ...
        param.ky_param, param.opt_param, param.iokr_param, param.debug_param.verbose);
    
    %--------------------------------------------------------------
    % Write out the IOKR Model
    %--------------------------------------------------------------
    outputDir = fullfile (outputDir, 'iokr_models');
    if ~ exist (outputDir, 'dir')
        mkdir (outputDir);
    end % if
    save (fullfile (outputDir, strcat (MP_IOKR_Defaults.param2str (param), '.mat')), ...
        'model', 'param');
end
#!/bin/bash

#--SBATCH --partition=debug --time=00:15:00
#--SBATCH --cpus-per-task=10 --mem-per-cpu=3000

#SBATCH --partition=batch --time=02:00:00
#SBATCH --cpus-per-task=20 --mem-per-cpu=1500
#SBATCH --nodes=1

module load matlab

command_str="\
    addpath (genpath ('/scratch/cs/kepaco/bache1/code/iokr-mp'));                                                                       \
    param = struct();                                                                                                                   \
    param.ky_param = struct ('representation', 'kernel', 'type', 'gaussian', 'base_kernel', 'linear', 'param_selection', 'entropy');    \
    param.iokr_param = struct ('mkl', 'alignf');                                                                                        \
    input_dir  = '/scratch/cs/kepaco/bache1/data/metabolite-identification/ACML17/';                                                    \
    output_dir = '/scratch/cs/kepaco/bache1/data/metabolite-identification/ACML17/results/iokr/nonlinear_output_kernel/';               \
    run_IOKR_ACML17 (input_dir, output_dir, param);                                                                                     \
    exit;"

matlab_multithread -nodisplay -r "$command_str"


    #param.debug_param = struct ('isDebugMode', true, 'n_debug_set', 250, 'randomSeed', 10);                                             \
    #param.ky_param = struct ('representation', 'feature', 'type', 'linear', 'base_kernel', 'linear', 'param_selection', 'cv');  \

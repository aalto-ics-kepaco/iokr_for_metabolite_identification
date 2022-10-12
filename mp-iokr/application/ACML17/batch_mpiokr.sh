#!/bin/bash

# -- SBATCH --partition=debug --time=00:15:00
# -- SBATCH --cpus-per-task=10 --mem-per-cpu=3000

#SBATCH --partition=batch --time=48:00:00
#SBATCH --cpus-per-task=20 --mem-per-cpu=10000
#SBATCH --nodes=1

module load matlab

command_str="\
    addpath (genpath ('/scratch/cs/kepaco/bache1/code/iokr-mp'));                                                                       \
    param = struct();                                                                                                                   \
    param.ky_param = struct ('representation', 'kernel', 'type', 'gaussian', 'base_kernel', 'tanimoto', 'param_selection', 'entropy');  \
    param.mp_iokr_param = struct ('mkl', 'unimkl');                                                                                     \
    param.data_param = struct();                                                                                                        \
    param.data_param.repetition      = ${1};                                                                                            \
    param.data_param.selection_param = struct ('strategy', 'random', 'perc', 1, 'inclExpCand', false);                                  \
    param.debug_param = struct ('randomSeed', ${SLURM_JOBID});                                                                  \
    input_dir  = '/scratch/cs/kepaco/bache1/data/metabolite-identification/ACML17/';                                                    \
    output_dir = '/scratch/cs/kepaco/bache1/data/metabolite-identification/ACML17/results/mp-iokr/nonlinear_output_kernel/';            \
    run_MP_IOKR_ACML17 (input_dir, output_dir, param);                                                                                  \
    exit;"

# echo $command_str
matlab_multithread -nodisplay -r "$command_str"


    # param.debug_param = struct ('isDebugMode', true, 'n_debug_set', 250, 'randomSeed', 10);                                             \

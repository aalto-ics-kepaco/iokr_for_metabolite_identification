# Preprocessing steps to train IOKR model for CSI-FingerID
Relevant sub-directories in the CSI-FingerID working directory:
- ```spectra```: MS/MS spectra / compounds used for training 
- ```kernels```: Kernel matrices containing base on all MS/MS spectra in 
                 ```spectra```
- ```fingerprints```: Fingerprint files for all compounds
- ```models_iokr```: IOKR model (as mat-file and binary) trained using *all* 
                     training data
- ```crossvalidation_models_iokr```: IOKR models (as mat-file and binary) trained
    for each pre-defined cross-validation fold (see ```crossvalidation_folds.txt```)

## The ```train_models_for_csifingerid``` script
This script can be used to train an IOKR model using *all* training data and get
the reclassficiation performance, i.e. predicted ranks and top-k accuarcy.

Main steps done in this script:

1. Load ```spectra/compound_info.mat```. This file contains the information about 
   each spectra / compound in the training set: 

| Column name | Description | 
| ---- | ---- |
| mol_id | basename of the spectrum-file, considered as compund id | 
| inchi_key_1 | first 14 characters (2D information) if the InChI-key | 
| inchi | 2D-InChI | 
| molecular_formula | molecular formula used to choose the candidate set | 
| fp_full | full fingerprint (type: sparse) | 
| fp_masked | fingerprint (type: sparse) masked using ```fingerprints/fingerprints.mask```| 

2. Load ```kernels/kernels.mat```
    * All kernels stored in ```kernels``` are used for training.
    * The kernels must have the file-extension '.txt', i.e. ```PPKr.txt```
    * If the kernels have not been yet converted into mat-files, it will be done
      automatically.
3. Train the model
    * Output-kernel parameters: Tanimoto-Gaussian Kernel, gamma selected using
      the entropy criteria
    * Input-kernels are combined using ALIGNF
4. Write out the model in binary format (```write_out_iokr_model_as_binary_files.m```)
    * All data is stored as double-precission binary (16byte, double) in little-ending
    * For each input kernel: 
        * Vector of column means: ```KX_mean_*.bin```
        * Kernel diagonal (centered): ```KX_diag_c_*.bin```
        * MKL-weight: ```mkl_w_*.bin```
    * Cholesky decomposition of __(lambda_opt * I + KX_train)__ in lower-triangular matrix
        * Vector contains *only* the elements of the lower-triangular (+ diagonal) 
          matrix.
        * The values are extracted row-wise.
        * ```C.bin```
    * For the output kernel:
        * Vector of column means: ```KY_mean.bin```
        * Kernel diagonal (centered): ```KY_diag_c.bin```
        * Gamma parameter for the Tanimoto-Gaussian kernel: ```KY_gamma.bin```


## Prepare the fingerprints
- The (unmasked) fingerprint corresponding to each spectra (compound) is stored 
  fingerprints.
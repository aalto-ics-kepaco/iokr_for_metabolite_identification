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

## Prepare the fingerprints
- The (unmasked) fingerprint corresponding to each spectra (compound) is stored 
  fingerprints.
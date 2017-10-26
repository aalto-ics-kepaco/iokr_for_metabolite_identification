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

1. Load the ```compound_info.mat``` (from ```spectra```)
    * This file contains the information about each spectra / compound in the 
      training set: 
| mol_id | basename of the spectrum-file, considered as kind of compund id | 


## Prepare the fingerprints
- The (unmasked) fingerprint corresponding to each spectra (compound) is stored 
  fingerprints.
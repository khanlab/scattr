# Output Files

After running the workflow, the `/path/to/output/dir` folder will contain a 
hidden `.logs` folder for troubleshooting, as well as additional folders 
associated with the tools that were used to generate files (e.g. `freesurfer`, 
`labelmerge`), but for most purposes, all the primary outputs of interest will 
be in the `mrtrix` directory with the following structure:

```
mrtrix/
└── sub-{subject}
    ├── dti
    ├── dwi
    ├── response
    └── tractography
```

Briefly, `dti` contains processed diffusion tensor images along with 
quantitative maps (e.g. `FA`, `MD`, etc.), `dwi` contains images that were
converted from Nifti (`.nii.gz`) to the MRtrix imaging format (`.mif`), 
`response` contains subject-specific response functions, and `tractography`
contains tractogram files.

## dti

This folder contains processed diffusion tensor images (DTI) that were used to 
compute quantitative maps (stored in the MRtrix imaging format - `.mif`), also 
found within the folder:

```
sub-{subject}
└── dti
    ├── sub-{subject}_model-dti_ad.mif
    ├── sub-{subject}_model-dti_fa.mif
    ├── sub-{subject}_model-dti_md.mif
    ├── sub-{subject}_model-dti_rd.mif
    └── sub-{subject}_model-dti_tensor.mif
```

As per the BIDS extension proposal, `model-dti` denotes to the DTI model used to
process the data, and the suffix (`ad`, `fa`, `md`, `rd`, `tensor`) describes 
the imaging data.

## dwi

This folder contains the input diffusion weighted image (dwi) and brain mask 
converted from Nifti (`.nii.gz`) to the MRtrix imaging format (`.mif`) used in
downstream processing. It also contains an intensity normalized dwi image used
for DTI:

```
sub-{subject}
└── dwi
    ├── sub-{subject}_brainmask.mif
    ├── sub-{subject}_desc-normalized_dwi.mif
    └── sub-{subject}_dwi.mif
```

## response

This folder contains the subject-specific response functions that were estimated
from the input diffusion data, as well as both the normalized and unnormalized
versions of the estimated fibre orientation distribution (FOD) maps. Both 
response functions and FOD maps are estimated for the three different tissue
types: white matter (wm), grey matter (gm), and corticospinal fluid (csf):

```
sub-{subject}
└── response
    ├── sub-{subject}_desc-csf_model-csd_fod.mif
    ├── sub-{subject}_desc-csf_model-csd_fodNormalized.mif
    ├── sub-{subject}_desc-csf_response.txt
    ├── sub-{subject}_desc-gm_model-csd_fod.mif
    ├── sub-{subject}_desc-gm_model-csd_fodNormalized.mif
    ├── sub-{subject}_desc-gm_response.txt
    ├── sub-{subject}_desc-wm_model-csd_fod.mif
    ├── sub-{subject}_desc-wm_model-csd_fodNormalized.mif
    └── sub-{subject}_desc-wm_response.txt
```

## tractography

This folder contains generated tractogram data and associated files describing
how streamlines connect different subcortical structures:

```
sub-{subject}
└── tractography
    ├── sub-{subject}_desc-filteredsubcortical_nodeAssignment.txt
    ├── sub-{subject}_desc-filteredsubcortical_nodeWeights.csv
    ├── sub-{subject}_desc-filteredsubcortical_tractography.tck
    ├── sub-{subject}_desc-iFOD2_muCoefficient
    ├── sub-{subject}_desc-iFOD2_tckWeights.txt
    ├── sub-{subject}_desc-iFOD2_tractography.tck
    ├── sub-{subject}_desc-subcortical_nodeAssignment.txt
    ├── sub-{subject}_desc-subcortical_nodeWeights.csv
```

The `desc-filteredsubcortical` entity pair describes associated tractogram 
files that have been filtered to only pass through WM and connect two GM 
structures, while the `desc-subcortical` entity pair denotes files associated
with the unfiltered subcortical connectome. Whole-brain tractogram associated
files are denoted by `desc-iFOD2`, which describes the tractography algorithm
used to perform tractograpy.

_Note: As the BIDS extension proposal has not been finalized, the naming of 
such files may be subject to change_

## Additional Files

The top-level `/path/to/output/dir` contains additional files / folders:
```
/path/to/output/dir
├── ...
├── config
├── .logs
├── .snakebids
└── .snakemake
```

The `config` folder, along with the hidden `.snakebids` and `.snakemake` folders 
contain a record of the code and parameters used, and paths to the inputs.

Workflow steps that write logs to file are stored in the hidden `.logs` 
subfolder, with the file names based on the tools used (e.g. `mrtrix`) and rule 
wildcards (e.g. `subject`). 

If the app is run in workflow mode (`--workflow-mode` / `-W`), which enables 
direct use of the Snakemake CLI to run scattr, output folders (e.g. `work`) will
be placed in a `results` folder.
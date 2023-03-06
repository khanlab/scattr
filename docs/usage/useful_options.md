# Running SCATTR on your data

This section goes over the command line options you will find most useful when
running SCATTR on your dataset, along with describing some issues you may face.

_Note: Please first refer to the simple example in the 
[Installation](https://scattr.readthedocs.io/en/stable/getting_started/installation.html) 
section, which goes over running SCATTR on a test dataset and the essential 
required options._

## Freesurfer License

To perform Freesurfer-related processing (e.g. thalamus segmentation), 
Freesurfer is directly invoked. As such, a Freesurfer license is required to 
perform such steps in the workflow. By default, SCATTR attempts to use the 
Freesurfer license saved in the environment variable `FS_LICENSE`.
Alternatively, the path to the Freesurfer license may be passed along by 
invoking the `--fs-license` parameter:

```
--fs-license /path/to/fs_license
```

## Including / excluding subjects to process
By default, SCATTR will run on **all** subjects in the dataset. If you wish to 
run on only a subset of subjects, you can use the `--participant-label` flag:

```
--participant-label 001
```

which would only run on `sub-001. You can add additional subjects by passing
a space-separated list to this option:

```
--participant-label 001 002
```

which would run for `sub-001` and `sub-002`.

Similarly, subjects can be excluded from processing using the 
`--exclude-participant-label` flag.

## Alternate Freesurfer / derived-diffusion data locations

By default, SCATTR attempts to locate Freesurfer and derived diffusion data 
locations. Users can overwrite these options, by passing along the actual 
location of each using either `--freesurfer_dir` or `--dwi_dir`, respectively.

```
--freesurfer_dir /path/to/fs_dir --dwi_dir /path/to/dwi_dir
```

## Pre-generated average response function
In some cases, an average response function may have already been separately 
generated, which is then used for downstream processing (e.g. generated from
controls and applied to a patient population). To use a pre-generated average
resposne function, the location of the directory containing the associated files
can be passed along using `--responsemean_dir`:

```
--responsemean_dir /path/to/average_response_dir
```

## Tractography on network storage

Performing tractography can require millions of reads and writes to the storage 
system in order to update the file on-the-fly. This process can be extremely 
slow on network storages. To help with this, SCATTR always reads and writes the
tractography to a temporary location (e.g. `/tmp`) before copying the output to
the final output, significantly improving the time it takes for tractography to 
be generated. On a network system, you may be unable to write to `/tmp`. An 
alternative on systems with `SLURM` workload managers is to invoke 
`--slurm_tmpdir`, which requests that the workflow write to the local temporary
storage system (e.g. `/localscratch`) instead of the network temporary storage.

## BIDS Parsing limitations

SCATTR uses Snakebids, which makes use of pybids to parse a [BIDS-compliant
dataset](https://bids.neuroimaging.io/). However, because of the way Snakebids
and Snakemake operate, one limitation is that the input files in your BIDS 
dataset needs to be consistent in terms of what optional BIDS entities exist in
them. We can use the acqusition (`acq`) entity as an example. SCATTR should have
no problem parsing the following dataset:

```
PATH_TO_BIDS_DIR/
└── dataset_description.json
└── sub-001/
    └── anat/
        ├── sub-001_acq-mprage_T1w.nii.gz
└── sub-002/
    └── anat/
        ├── sub-002_acq-spgr_T1w.nii.gz
...
```

as the path (with wildcards) will be interpreted as 
`sub-{subject}_acq-{acq}_T1w.nii.gz`.

However, the following dataset will raise an error:

```
PATH_TO_BIDS_DIR/
└── dataset_description.json
└── sub-001/
    └── anat/
        ├── sub-001_acq-mprage_T1w.nii.gz
└── sub-002/
    └── anat/
        ├── sub-002_T1w.nii.gz
...
```

because two distinct paths (with wildcards) would be found for T1w images:
`sub-{subject}_acq-{acq}_T1w.nii.gz` and `sub-{subject}_T1w.nii.gz`.

Similarly, you could not have some subjects with the `ses` identifier, and some
subjects without it. There will soon be added functionality in Snakebids to 
filter out extra files, but for now, if your dataset has these issues, you will
need to rename or remove extraneous files.

More example of possible BIDS-compliant datasets can be found in 
`scattr/test/data`.
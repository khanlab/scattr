# Frequently Asked Questions (FAQs)

### 1. Why do I get the error `No input images found for <modality>`?

The workflow was unable to find the required input files from the input BIDS
directory. This can happen if:

* Singularity or Docker cannot access your input directory. For Singularity,
ensure your 
[Singularity options](https://docs.sylabs.io/guides/3.1/user-guide/cli/singularity_run.html) 
are appropriate, in particular `SINGULARITY_BINDPATH`. For Docker, ensure you 
are mounting the correct directory with the `-v` flag described in the 
[Getting Started](https://scattr.readthedocs.io/en/stable/getting_started/docker.html)
section.
* SCATTR does not recognize your BIDS-formatted input images. This can occur if,
for example, T1w images are labelled with the suffix `t1w.nii.gz` instead of 
`T1w.nii.gz` as per [BIDS specifications](https://bids.neuroimaging.io/specification.html).
SCATTR makes use of [pybids](https://github.com/bids-standard/pybids) to parse
the dataset, so we suggest using the [BIDS Validator](https://bids-standard.github.io/bids-validator/) to ensure your dataset
has no errors.

If you passed Freesurfer or diffusion derivative locations (`--freesurfer-dir` 
or `--dwi-dir`), you may get get this message as the files are a different
location than in the input directory provided.
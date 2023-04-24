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
the dataset, so we suggest using the 
[BIDS Validator](https://bids-standard.github.io/bids-validator/) to ensure 
your dataset has no errors.

If you passed Freesurfer or diffusion derivative locations (`--freesurfer-dir` 
or `--dwi-dir`), you may get get this message as the files are a different
location than in the input directory provided.

### 2. What if I want to use merge two different atlases?
<!-- Update with readthedocs link when available -->
`SCATTR` uses [labelmerge](https://github.com/khanlab/labelmerge) to combine
labels from different atlases. To that end, we have exposed the `labelmerge` 
command-line arguments, enabling substitution of different atlases. These
arguments are pre-pended with `--labelmerge_<labelmerge_arg>`. We encourage you
to check out the `labelmerge` documention (linked above), as well as `SCATTR`'s 
help command to see all available options. 

*NOTE:* If using a different atlas, you will need to ensure that the associated
metadata is available (according to the 
[BIDS specification](https://bids-specification.readthedocs.io/en/stable/05-derivatives/03-imaging.html#segmentations)) 
for each subject part of the input dataset and in the same space as the T1w 
image.

### 3. I only want to use labels from a single atlas.
If you only want to use a labels from a single atlas (e.g. skip labelmerge), you
can invoke the `--skip_labelmerge` argument. If you are using this flag, you can
combine this with the `--labelmerge_base_dir` and `--labelmerge_base_desc` 
arguments to provide your own atlas. Similar to using custom atlases (from _2_),
it is expected that the imaging and associated metadata is available for all 
subjects in the dataset and that all files follow the BIDS specification.

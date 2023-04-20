# Running SCATTR with Singularity

## Pre-requisites
1. Singularity / Apptainer is is installed on your system. For more info, see
the detailed [Apptainer install instructions](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages).
1. The following command-line tools are installed:
    * wget
1. Sufficient disk-space (rough estimate)
    * in your `/tmp` folder (>30GB) to build the container
    * in your working folder to store the container (~20GB)
    * for SCATTR outputs (~40GB per subject using default parameters)
1. Sufficient CPU and memory - the more you have, the faster it will run and 
the more streamlines that can be estimated. We recommand at least 8 CPU cores 
and 64GB memory if using default parameters.

## First time setup
Pull the container. This can be done from DockerHub, but requires a large 
amount of disk space in your `/tmp` folder, since it has to convert from a 
Docker container to a Singularity/Apptainer container. The example below pulls
the latest versioned container (replace `latest` with `vX.X.X` for a specific
version).

```
singularity pull docker://khanlab/scattr:latest
```
_Note: If you encounter any errors pulling the container from DockerHub, it may
be because you are running out of disk space in your cache folders. You can 
change these locations by setting environment variables, however, using a 
network file system for the folders may result in poor performance:_
```
export SINGULARITY_CACHEDIR=/YOURDIR/.cache/singularity
```


Run SCATTR without any arguments to print the short help:

```
singularity run -e khanlab_scattr_latest.sif
```

Use the `-h` option to get a detailed help listing:

```
singularity run -e khanlab_scattr_latest.sif -h
```

Note that all the Snakemake command-line options are also available in SCATTR,
and can be listed with `--help-snakemake`:

```
singularity run -e khanlab_scattr_latest.sif --help-snakemake
```

## Running an example

We will use the `test` folder found from the 
[Github repository](https://github.com/khanlab/scattr/tree/main/test/) to
demonstrate an example of how to run SCATTR:

```
singularity run -e khanlab_scattr_latest.sif test/data/bids test/data/derivatives participant --fs-license test/fs_license --force-output -n
```

### Explanation

Everything prior to the container (`khanlab_scattr_latest.sif`) are arguments
to Singluarity / Apptainer, and after are to SCATTR itself. The first three arguments to 
SCATTR (as with any BIDS App) are the input folder (`test/data/bids`), the 
output folder (`test/data/derivatives`), and the analysis level (`participant`).
The `participant` analysis level is used in SCATTR to perform further 
participant-level processing of Freesurfer (thalamus segmentation), external 
atlases (combining segmentations) and diffusion derived data (estimation of 
fibre orientation distributions). This includes estimating an average response 
function from input data. The `--fs-license` argument allows for specification
of the location of the required Freesurfer license file if not already 
specified in the `FS_LICENSE` environment variable. Note, that this is
required to perform any Freesurfer-related processing. The `--force-output` 
flag is a Snakemake argument that is invoked to allow for writing of output file
to already existing folders - in this case, for thalamus segmentations via 
Freesurfer. We also used the `--dry-run/-n` option to print out what would run,
without actually running the workflow.

When you run the above command, a long listing will print out, describing all 
the rules that will be run. We can also have a shell command used for each rule
printed to screen using the `-p` Snakemake option

```
singularity run -e khanlab_scattr_latest.sif test/data/bids test/data/derivatives participant --fs-license test/fs_license --force-output -np
```

Now to actually run the workflow, we need to specify how many cores to use and 
leave out the dry-run option. The Snakemake `--cores` option tells SCATTR how
many cores to use. Using `--cores 8` means that SCATTR will only make use of 8 
cores at most. Generally speaking, you should use `--cores all`, so it can make 
maximal use of all available CPU cores it has access to on your system. This is 
especially useful if you are running multiple subjects.

Running the follow command (SCATTR on a single subject) may take up to ~36 hours
with 8 cores and default parameters, but could be much longer (several days) if 
you only have a single core.

```
singularity run -e khanlab_scattr_latest.sif test/data/bids test/data/derivatives participant --fs-license test/fs_license --force-output -p --cores all
```

_Note that you may need to adjust your 
[Singularity / Apptainer options](https://sylabs.io/guides/3.1/user-guide/cli/singularity_run.html) 
to ensure the container can read and write to your input and output directories, 
respectively. You can bind paths easily by setting an environment variable, 
e.g. if you have a `/project` folder that contains your data, you can add it to
the `SINGULARITY_BINDPATH` so it is available when you are running a container:_

```
export SINGULARITY_BINDPATH=/data:/data
```

After this completes, you have additional folders in your output folder,
`test/data/derivatives`, for the one subject.

### Exploring different options

If you alternatively want to run SCATTR using a pre-defined average response 
function, you can use the `--responsemean_dir` flag to specify the location to
where the average response function is located. 

```
singularity run khanlab_scattr_latest.sif test/data/bids test/data/derivatives/ participant --responsemean_dir test/data/derivatives/mrtrix/avg --fs-license test/.fs_license -np --force-output
```

Other parameters exist, which may help to improve processing times at the 
expense of sensitivity / specificity (e.g. reducing the number of streamlines 
generated).
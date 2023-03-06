# Running SCATTR with Docker on Windows

_Note, these instructions you have Docker installed already on a Windows system.
Docker can also be run on Linux or MacOS with similar commands, but here, we 
will assume the default Windows CLI is being used.

## First time setup

Open your Windows Command Prompt by clicking the `Windows` button and typing
`cmd` and pressing the `Enter` on your keyboard. This is where you will enter 
your SCATTR commands. Feel free to make a new directory with `mkdir` or move to
a directory you would like to work out of with `cd. For this example, we will
work from:

```
cd c:\Users\username\Downloads
```

Pull the container (this will take some time and storage stage, but like an 
installation, it only needs to be done once and can be then be run on many 
datasets). The example below pulls the latest versioned container (replace 
`latest` with `vX.X.X` for a specific version).

```
docker pull khanlab/scattr:latest
```

Run SCATTR without any arguments to print hte short help:

```
docker run -it --rm khanlab/scattr:latest
```

Use the `-h` option to get a detailed help listing:

```
docker run -it --rm khanlab_scattr_latest.sif -h
```

_Note that all the Snakemake command-line options are also available in SCATTR,
and can be listed with `--help-snakemake`:

```
docker run -it --rm khanlab_scattr_latest.sif --help-snakemake
```

## Running an example

We will use the `test` folder found from the 
[Github repository](https://github.com/khanlab/scattr/tree/main/test/) via
`git clone` to the previously mentioned folder to demonstrate an example of 
how to run SCATTR

```
docker run -it --rm -v c:\Users\username\Downloads\scattr\test:\test khanlab_scattr_latest.sif /test/data/bids /test/data/derivatives participant --fs-license /test/fs_license --force-output -n
```

### Explanation

Everything prior to the container (`khanlab_scattr_latest.sif`) are arguments
to Docker and after are to SCATTR itself. The first three arguments to Docker
are to enable interactive mode (`-it`), run and subsequently remove the Docker
container upon completion (`--rm`) and mount the the directoty 
(`-v c:\Users\username\Downloads\scattr\test`) to a directory within the
container named `\test`. These are not specific to SCATTR, but are general ways
to use Docker. You may want to familiarize yourself with 
[Docker options](https://docs.docker.com/engine/reference/run/).

The first three arguments to SCATTR (as with any BIDS App) are the input folder 
(`test/data/bids`), the output folder (`test/data/derivatives`), and the 
analysis level (`participant`). The `participant` analysis level is used in 
SCATTR to perform further participant-level processing of Freesurfer (thalamus 
segmentation), external atlases (combining segmentations) and diffusion derived
data (estimation of fibre orientation distributions). This includes estimating 
an average response function from input data. The `--fs-license` argument allows
for specification of the location of the required Freesurfer license file if not 
already specified in the `FS_LICENSE` environment variable. Note, that 
this is required to perform any Freesurfer-related processing. The 
`--force-output` flag is a Snakemake argument that is invoked to allow for 
writing of output file to already existing folders - in this case, for thalamus 
segmentations via Freesurfer. We also used the `--dry-run/-n` option to print 
out what would run, without actually running the workflow.

When you run the above command, a long listing will print out, describing all 
the rules that will be run. We can also have a shell command used for each rule
printed to screen using the `-p` Snakemake option

```
docker run -it --rm -v c:\Users\username\Downloads\scattr\test:\test  khanlab_scattr_latest.sif /test/data/bids /test/data/derivatives participant --fs-license /test/fs_license --force-output -np
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
docker run -it --rm -v c:\Users\username\Downloads\scattr\test:\test /test/data/bids /test/data/derivatives participant --fs-license /test/fs_license --force-output -p --cores all
```

After this completes, you have additional folders in your output folder,
`c:\Users\username\Downloads\scattr\test\data\derivatives`, for the one subject.

### Exploring different options

If you alternative want to run SCATTR using a pre-defined average response 
function, you can use the `--responsemean_dir` flag to specify the location to
where the average response function is located. 

```
docker run -it --rm -v c:\Users\username\Downloads\scattr\test:\test /test/data/bids /test/data/derivatives/ participant --responsemean_dir /test/data/derivatives/mrtrix/avg --fs-license /test/.fs_license -np --force-output
```

Other parameters exist, which may help to improve processing times at the 
expense of sensitivity / specifity (e.g. reducing the number of streamlines 
generated).
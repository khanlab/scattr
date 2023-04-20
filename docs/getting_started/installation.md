# Installation
BIDS App for Structural Connectivity Applied To Targeted Regions (SCATTR)

## Requirements
* Docker (Intel Mac/Windows/Linux) or Singularity (Linux)
* For those wishing to contribute or modify the code, `pip install` or `poetry
install` are also available (Linux), but will still require Singularity to 
handle some dependencies. See 
[Contributing to SCATTR](https://scattr.readthedocs.io/en/stable/contributing/contributing.html).
* _Note: Apple ARM-based chips (e.g. M1, M2, etc.) are **not currently
supported**. We do not have a Docker arm64 container yet._

### Notes
* Inputs to SCATTR should be a BIDS dataset, including processed Freesurfer
and DWI derivatives (which can be stored separately). 

## Docker on Windows / Mac (Intel) / Linux

The SCATTR BIDS App is available on DuckerHub as versioned releases.
Instructions can be found in the [Docker](https://scattr.readthedocs.io/en/stable/getting_started/docker.html) documentation page.

### Pros
* Compatible with non-Linux systems
* All dependencies are in a single container

### Cons
* Typically not possible on shared machines
* Cannot use Snakemake cluster execution profiles
* Cannot edit code

## Singularity Container

The same Docker container can also be used with Singularity (now Apptainer).
Instructions can be found in the [Singularity / Apptainer](https://scattr.readthedocs.io/en/stable/getting_started/singularity.html) documentation page.

### Pros
* All dependencies are in a single container, stored as a single file (.sif)
* Compatible on shared systems with Singularity installed

### Cons
* Cannot use Snakemake cluster execution profiles
* Cannot edit code

## Python Environment with Singularity Dependencies

Instructions can be found in the [Contributing](https://scattr.readthedocs.io/en/stable/contributing/contributing.html) documentation page.

### Pros
* Flexibility to modify code

### Cons
* Only compatible on systems with Singularity for external dependencies
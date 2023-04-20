# Workflow Details

This section describes the SCATTR workflow (i.e. steps taken to produce 
intermediate and final files). SCATTR is a 
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow, and thus a 
directed acyclic graph (DAG) that is automatically configured based on a set of 
rules.

## Overall workflow

Below is an image exhibiting the workflow DAG. Each rounded rectangle in 
the DAG represents a rule (i.e. some code or script that produces an output), 
with arrows representing the connections (i.e. inputs / outputs) to these rules.

<img src="https://raw.githubusercontent.com/khanlab/scattr/main/docs/workflow/dag.png" width="800px">

Although it may look complex, it is also organized into groups of rules, each
representing the different phases of the workflow. Each grouped set of rules 
also exist in separate rule files, which can be found under the 
[rules sub-directories](https://github.com/khanlab/scattr/tree/main/scattr/workflow/rules) 
in the workflow source. For example, the [`freesurfer.smk`](https://github.com/khanlab/scattr/tree/main/scattr/workflow/rules/freesurfer.smk)
file contains rules associated with further processing using Freesurfer 
(e.g. thalamus segmentation), and these are grouped together in the above 
diagram by a blue rectangle labeled `freesurfer`.

The main phases of the workflow are described in the sections below, zooming in 
on the rules within each blue rectangle.

### Freesurfer

Segmentation of the thalamus into its subnuclei is performed via the Freesurfer
script. As the output from this segmentation is in the Freesurfer file format 
(.mgz), a conversion to nifti (.nii.gz) is performed. Once converted, a 
transformation is computed to align the segmentations from the Freesurfer space 
to the subject's native space. The computed transformation is applied to the 
segmentations and used for further downstream processing.

![Freesurfer workflow](https://raw.githubusercontent.com/khanlab/scattr/main/docs/workflow/freesurfer_dag.png)

### Zona BB Subcortex

SCATTR also uses [labelmerge](https://github.com/khanlab/labelmerge) to combine
segmentations from varying sources. This allows for the previously segmented
thalamic nuclei to be combined with an atlas of other subcortical structures. 
This newly combined atlas is used for two separate purposes:

1. To individual binarized masks for the structures of interest. This is used to
identify connections terminating within these structures.
1. To create a convex hull (including the brainstem) allow for identification
of only those connections within the larger region of interest (e.g. 
subcortical region around the subcortical structures)

![Subcortex workflow](https://raw.githubusercontent.com/khanlab/scattr/main/docs/workflow/subcortex_dag.png)

### Mrtpipelines

[MRtrix3](https://www.mrtrix.org) is used in the workflow for 
tractography-related processing. These rules first convert the necessary nifti
images to the MRtrix file format, `.mif`. Similar to the subcortex workflow, two
different processes are performed:

1. Diffusion tensor imaging is performed to calculate quantitative 
maps of fractional anisotropy, mean diffusivity, radial diffusvity and 
axial diffusivity.
1. Individual subject response functions for the three tissue types, white 
matter, grey matter, and corticospinal fluid, are first estimated, and an 
average response function is then computed. Alternatively, a pre-computed 
response function can be supplied to be used in the downstream rules. From the 
average response functions, fibre orientation distribution (FOD) maps are 
computed for individual subjects. FOD maps are normalized to correct for 
the effects of (residual) intensity inhomogeneties. Once normalized, whole-brain
tractogaphy is computed from the FOD maps, before spherical-deconvolutional 
filtering of tractograms (`tcksift2`) is performed to weight the tractogram to
match the underlying diffusion signal. Connections between specific structures
of interest are identified, and further filtered to retain only those passing
through the white matter (`tck2connectome`, `connectome2tck`, 
`filter_combine_tck`, and `filtered_tck2connectome`).

![Mrtrix workflow](https://raw.githubusercontent.com/khanlab/scattr/main/docs/workflow/mrtpipelines_dag.png)
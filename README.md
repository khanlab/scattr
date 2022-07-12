## Deep Brain Structural Connectivity (DBSC)

**Development in-progress** (conversion from shell scripts)

DBSC is a workflow that maps connections between subcortical structures.
It uses existing tools like `prepdwi`, `Freesurfer`, and `Mrtrix3`, as well
makes use of existing parcellations.

This workflow was used to process the data for the analysis from 
[`hcp_subcortical_repo`](https://github.com/kaitj/hcp_subcortical_repro).

### Notes

* Original workflow had worked on HCP data which had transforms readily available
Transforms need to be computed to standard spaces used (MNI2009b is sufficient)
* Current workflow assumes `prepdwi` and `Freesurfer` has already been run.
Future implementation will include option to run these workflows.
* The SnakeBIDS workflow has not yet been tested.
Original workflow was run using scripts similar to those in `resources/example_scripts`

### Relevant Papers

Kai, J, Khan AR, Haast RAM, Lau JC. Mapping the subcortical connectome using in vivo diffusion MRI: feasibility and reliability. bioRxiv 2022.03.28.485689. doi: [10.1101/2022.03.28.485689](https://doi.org/10.1101/2022.03.28.485689)

## Deep Brain Structural Connectivity (DBSC)

**Development in-progress** (conversion from shell scripts)

DBSC is a BIDS App that performs tractography to identify 
connections between subcortical structures (i.e. map the subcortical 
connectome), making use of existing neuroimaging tools like `prepdwi`, 
`Freesurfer`, and `Mrtrix3`.

This workflow was used to process the data for the analysis from 
[`hcp_subcortical_repo`](https://github.com/kaitj/hcp_subcortical_repro).

### Contributing
Clone the git repository. DBSC dependencies are managed with Poetry 
(version 1.2.x), which you'll need installed on your machine. 
You can find instructions on the Poetry 
[website](https://python-poetry.org/docs/). 

Then, setup the development environment with the following commands:

```
poetry install
poetry run poe setup
```

DBSC uses poethepoet as a task runner. 
You can see what commands are available by running:

```
poetry run poe
```

If you wish, you can also run poe [command] directly by installing poethepoet 
on your system. Follow the install instructions at the link above.

DBSC uses pre-commit hooks (installed via the poe setup command above) to lint 
and format code (we use black, isort, flake8). By default, these hooks are 
run on every commit. Additionally, please run the following task:

```
poetry run poe quality
```

Please be sure they all pass before making a PR.

### Notes

* Original workflow had worked on HCP data which had transforms readily 
available. Transforms need to be computed to standard spaces used 
(MNI2009b is sufficient).
* Current workflow assumes `prepdwi` and `Freesurfer` has already been run.
Future implementation will include option to run these workflows.
* Original workflow was run using scripts similar to those in 
`resources/example_scripts`

### Relevant Papers

* Kai, J., Khan, A.R., Haast, R.A.M., Lau, J.C. (2022). 
Mapping the subcortical connectome using in vivo diffusion MRI: feasibility 
and reliability. Terra incognita: diving into the human subcortex, 
special issue of NeuroImage. 
doi: [10.1016/j.neuroimage.2022.119553](https://doi.org/10.1016/j.neuroimage.2022.119553).
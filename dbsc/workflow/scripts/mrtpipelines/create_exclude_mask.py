#!/usr/bin/env python
from pathlib import Path

from snakebids import bids
from snakemake.io import expand
from snakemake.shell import shell

def create_exclude_mask(
    container,
    base_dir,
    subcortical_seg, 
    num_labels, 
    out_dir, 
    wildcards,
    subj_wildcards,
    threads,
):
    # Create output directory
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    lZI = bids(
        root=base_dir,
        datatype="anat",
        desc="21",
        suffix="mask.mif",
    )
    
    rZI = bids(
        root=base_dir,
        datatype="anat",
        desc="21",
        suffix="mask.mif",
    )

    # Build exclude mask list
    exclude_mask_list = []
    for node1, node2 in [*zip(input.roi1, input.roi2)]:
        exclude_mask_list.extend(
            expand(
                bids(
                    root=base_dir,
                    datatype="exclude_mask",
                    desc=f"from{node1}-{node2}",
                    suffix="mask.mif",
                    **subj_wildcards,
                ),
                **wildcards,
            ),
        )

    # Create masks - parallelizes within a single job
    shell("singularity run {container} parallel --jobs {threads} mrcalc {subcortical_seg} 0 -neq {{1}} -sub {{2}} -sub {{lZI}} -sub {{rZI}} {{3}}} ::: {roi1} :::+ {roi2} :::+ {exclude_mask_list}")

    shell("rm -r {mask_dir}")

if __name__ == "__main__":
    create_roi_mask(
        container=snakemake.params.container, #noqa: F821
        base_dir=snakemake.params.base_dir, #noqa: F821
        roi1=snakemake.input.roi1, #noqa: F821
        roi2=snakemake.input.roi2, #noqa: F821
        subcortical_seg=snakemake.input.subcortical_seg, #noqa: F821
        mask_dir=snakemake.params.mask_dir, #noqa: F821
        out_dir=snakemake.output.out_dir, #noqa: F821
        wildcards=snakemake.wildcards, #noqa: F821
        subj_wildcards=snakemake.params.subj_wildcards, #noqa: F821
        threads=snakemake.threads, #noqa: F821
    )

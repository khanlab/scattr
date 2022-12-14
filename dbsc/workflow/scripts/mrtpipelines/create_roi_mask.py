#!/usr/bin/env python
from pathlib import Path

from snakebids import bids
from snakemake.io import expand
from snakemake.shell import shell

def create_roi_mask(
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

    # Read number of labels
    with open(num_labels, "r") as f:
        num_labels = int(f.read().strip())

    # Build mask list
    roi_mask_list, label_list = [], []
    for node_idx in range(1, num_labels+1):
        label_list.append(node_idx)
        roi_mask_list.extend(
            expand(
                bids(
                    root=base_dir,
                    datatype="roi_masks",
                    desc=node_idx, 
                    suffix="mask.mif",
                    **subj_wildcards,
                ), 
                **wildcards,
            )
        )

    out = ' '.join(mask for mask in roi_mask_list)

    # Get individual masks - parallelize within job
    shell("singularity run {container} parallel --jobs {threads} mrcalc {subcortical_seg} {{1}} -eq {{2}} ::: {label_list} :::+ {out}")


if __name__ == "__main__":
    create_roi_mask(
        container=snakemake.params.container, #noqa: F821
        base_dir=snakemake.params.base_dir, #noqa: F821
        subcortical_seg=snakemake.input.subcortical_seg, #noqa: F821
        num_labels=snakemake.input.num_labels, #noqa: F821
        out_dir=snakemake.output.out_dir, #noqa: F821
        wildcards=snakemake.wildcards, #noqa: F821
        subj_wildcards=snakemake.params.subj_wildcards, #noqa: F821
        threads=snakemake.threads, #noqa: F821
    )

#!/usr/bin/env python
from pathlib import Path

from snakebids import bids
from snakemake.io import expand
from snakemake.shell import shell

import numpy as np

def create_exclude_mask(
    container,
    base_dir,
    roi1, 
    roi2,
    subcortical_seg,
    num_labels,
    mask_dir, 
    out_dir, 
    wildcards,
    threads,
):
    # Create output directory
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Get number of labels
    with open(num_labels, "r") as f:
        num_labels = int(f.read().strip())

    node_pairs = np.triu_indices(num_labels, k=1)

    # Zona labels (currently misisng subjects)
    lZI=expand(
        bids(
            root=base_dir,
            datatype="roi_masks",
            desc="21",
            suffix="mask.mif",
            **wildcards,
        )
    )  
    rZI=bids(
        root=base_dir,
        datatype="roi_masks",
        desc="22",
        suffix="mask.mif",
        **wildcards,
    )

    # Create mask sequentially
    out_mask = expand(
        bids(
            root=base_dir,
            datatype="exclude_mask",
            desc="{desc}",
            suffix="mask.mif",
            **wildcards,
        ),
        desc=[
            f'from{node_pairs[0][idx]+1}-{node_pairs[1][idx]+1}' for idx in range(
                len(node_pairs[0])
            )
        ],
    )

    # Create masks - parallelizes within a single jobs
    shell("singularity run {container} parallel --citation --jobs {threads} -k mrcalc {subcortical_seg} 0 -neq {{1}} -sub {{2}} -sub {lZI} -sub {rZI} -sub {{3}} ::: {roi1} :::+ {roi2} :::+ {out_mask}")

if __name__ == "__main__":
    create_exclude_mask(
        container=snakemake.params.container, #noqa: F821
        base_dir=snakemake.params.base_dir, #noqa: F821
        roi1=snakemake.input.roi1, #noqa: F821
        roi2=snakemake.input.roi2, #noqa: F821
        subcortical_seg=snakemake.input.subcortical_seg, #noqa: F821
        num_labels=snakemake.input.num_labels, #noqa: F821
        mask_dir=snakemake.params.mask_dir, #noqa: F821
        out_dir=snakemake.output.out_dir, #noqa: F821
        wildcards=snakemake.wildcards, #noqa: F821
        threads=snakemake.threads, #noqa: F821
    )
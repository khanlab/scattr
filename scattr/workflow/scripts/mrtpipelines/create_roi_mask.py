#!/usr/bin/env python
from pathlib import Path

from snakebids import bids
from snakemake.io import expand  # (used in commented command)
from snakemake.shell import shell


def create_roi_mask(
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
    with open(num_labels) as f:
        num_labels = int(f.read().strip())

    # Create masks sequentially (error with trying to do it in parallel)
    for idx in range(1, num_labels + 1):
        roi = bids(  # (used in shell call)
            root=base_dir,
            datatype="roi_masks",
            desc=f"{idx}",
            suffix="mask.mif",
            **wildcards,
        )

        shell(
            "mrcalc -nthreads {threads} "
            "{subcortical_seg} {idx} -eq {roi} -force"
        )

    # Build mask list and create in parallel
    # label_list = [idx for idx in range(1, num_labels+1)]
    # roi_mask_list = expand(
    #     bids(
    #         root=base_dir,
    #         datatype="roi_masks",
    #         desc="{desc}",
    #         suffix="mask.mif",
    #         **wildcards,
    #     ),
    #     desc=[idx for idx in range(1, num_labels+1)]
    # )
    # out = ' '.join(mask for mask in roi_mask_list)

    # shell(
    #     "singularity run {container} parallel --jobs {threads} -k "
    #     "mrcalc {subcortical_seg} {{1}} -eq {{2}} ::: {label_list} "
    #     ":::+ {out}"
    # )


if __name__ == "__main__":
    create_roi_mask(
        base_dir=snakemake.params.base_dir,
        subcortical_seg=snakemake.input.subcortical_seg,
        num_labels=snakemake.input.num_labels,
        out_dir=snakemake.output.out_dir,
        wildcards=snakemake.wildcards,
        subj_wildcards=snakemake.params.subj_wildcards,
        threads=snakemake.threads,
    )

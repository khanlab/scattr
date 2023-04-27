def aggregate_tck_files(wildcards):
    """Grab all files associated with unfiltered tck"""
    # Build list of files
    unfiltered_weights = expand(
        bids_tractography_out(
            datatype="unfiltered",
            desc="subcortical",
            suffix="tckWeights{suffix}.csv",
        ),
        **wildcards,
        allow_missing=True,
    )[0]

    suffix = glob_wildcards(unfiltered_weights).suffix

    unfiltered_tck = expand(
        bids_tractography_out(
            datatype="unfiltered",
            desc="subcortical",
            suffix="from{suffix}.tck",
        ),
        suffix=suffix,
        allow_missing=True,
    )

    return {
        "weights": expand(
            bids_tractography_out(
                datatype="unfiltered",
                desc="subcortical",
                suffix="tckWeights{suffix}.csv",
            ),
            suffix=suffix,
            allow_missing=True,
        ),
        "tck": unfiltered_tck,
    }


def aggregate_exclude_masks(wildcards):
    """Grab all exclude masks"""
    # Build list of files
    exclude_mask = expand(
        bids_anat_out(
            datatype="exclude_mask",
            desc="{desc}",
            suffix="mask.mif",
        ),
        **wildcards,
        allow_missing=True,
    )

    # Get desc wildcard
    desc = glob_wildcards(exclude_mask[0]).desc

    return expand(
        exclude_mask,
        desc=desc,
    )


def get_desc(wildcards):
    """Build param list"""
    # Get desc wildcards
    exclude_mask = expand(
        bids_anat_out(
            datatype="exclude_mask",
            desc="{desc}",
            suffix="mask.mif",
        ),
        **wildcards,
        allow_missing=True,
    )[0]
    desc = glob_wildcards(exclude_mask).desc

    return desc


def get_filtered_tck(wildcards):
    desc = get_desc(wildcards)
    # Create params lists
    return expand(
        bids_tractography_out(
            desc="{desc}",
            suffix="tractography.tck",
        ),
        **wildcards,
        desc=desc,
        allow_missing=True,
    )


def get_filtered_weights(wildcards):
    desc = get_desc(wildcards)
    return expand(
        bids_tractography_out(
            desc="{desc}",
            suffix="weights.csv",
        ),
        **wildcards,
        desc=desc,
        allow_missing=True,
    )


rule filter_combine_tck:
    input:
        unpack(aggregate_tck_files),
        rules.connectome2tck.output.output_dir,
        rules.create_exclude_mask.output.out_dir,
        filter_mask=aggregate_exclude_masks,
    params:
        filtered_tck=get_filtered_tck,
        filtered_weights=get_filtered_weights,
        filtered_tck_exists=bids_tractography_out(
            desc="from*", suffix="tractography.tck"
        ),
        filtered_weights_exists=bids_tractography_out(
            desc="from*", suffix="weights.csv"
        ),
        exclude_mask_dir=bids_anat_out(
            datatype="exclude_mask",
        ),
        unfiltered_tck_dir=str(
            Path(
                bids_tractography_out(
                    datatype="unfiltered",
                )
            ).parent
        ),
    output:
        combined_tck=bids_tractography_out(
            desc="filteredsubcortical",
            suffix="tractography.tck",
        ),
        combined_weights=bids_tractography_out(
            desc="filteredsubcortical",
            suffix="tckWeights.txt",
        ),
    threads: 1
    resources:
        tmp_dir=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            **wildcards,
        ),
        tmp_combined_tck=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            desc="filteredsubcortical",
            suffix="tractography.tck",
            **wildcards,
        ),
        tmp_combined_weights=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            desc="filteredsubcortical",
            suffix="tckWeights.txt",
            **wildcards,
        ),
        mem_mb=128000,
        time=60,
    log:
        bids_log(suffix="combineFiltered.log"),
    group:
        "tractography_update"
    container:
        config["singularity"]["scattr"]
    shell:
        """
        mkdir -p {resources.tmp_dir}

        parallel --jobs {threads} -k tckedit -exclude {{1}} \\
            -tck_weights_in {{2}} -tck_weights_out {{3}} {{4}} \\
            {{5}} ::: {input.filter_mask} :::+ {input.weights} :::+ \\
            {params.filtered_weights} :::+ {input.tck} :::+ \\
            {params.filtered_tck} || true

        tckedit {params.filtered_tck_exists} \\
            {resources.tmp_combined_tck} &> {log}

        cat {params.filtered_weights_exists} >> \\
            {resources.tmp_combined_weights}

        rsync -v {resources.tmp_combined_tck} \\
            {output.combined_tck} >> {log} 2>&1

        rsync -v {resources.tmp_combined_weights} \\
            {output.combined_weights} >> {log} 2>&1

        rm {params.filtered_tck_exists} {params.filtered_weights_exists}
        """


rule filtered_tck2connectome:
    input:
        weights=rules.filter_combine_tck.output.combined_weights,
        tck=rules.filter_combine_tck.output.combined_tck,
        subcortical_seg=rules.tckgen.input.subcortical_seg,
    params:
        radius=config["radial_search"],
    output:
        sl_assignment=bids_tractography_out(
            desc="filteredsubcortical",
            suffix="nodeAssignment.txt",
        ),
        node_weights=bids_tractography_out(
            desc="filteredsubcortical",
            suffix="nodeWeights.csv",
        ),
    threads: 32
    resources:
        tmp_dir=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            **wildcards,
        ),
        tmp_sl_assignment=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            desc="filteredsubcortical",
            suffix="nodeAssignment.txt",
            **wildcards,
        ),
        tmp_node_weights=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            desc="filteredsubcortical",
            suffix="nodeWeights.csv",
            **wildcards,
        ),
        mem_mb=128000,
        time=60 * 3,
    log:
        bids_log(desc="filtered", suffix="tck2connectome"),
    group:
        "tractography_update"
    container:
        config["singularity"]["scattr"]
    shell:
        """mkdir -p {resources.tmp_dir}

        tck2connectome -nthreads {threads} -zero_diagonal -stat_edge sum \\
            -assignment_radial_search {params.radius} \\
            -tck_weights_in {input.weights} \\
            -out_assignments {resources.tmp_sl_assignment} \\
            -symmetric {input.tck} {input.subcortical_seg} \\
            {resources.tmp_node_weights} &> {log}

        rsync {resources.tmp_sl_assignment} \\
            {output.sl_assignment} >> {log} 2>&1

        rsync {resources.tmp_node_weights} \\
            {output.node_weights} >> {log} 2>&1
        """

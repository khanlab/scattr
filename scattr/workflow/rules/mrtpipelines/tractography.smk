rule tckgen:
    """
    Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic 
    streamlines tractography by 2nd order integration over fibre orientation 
    distributions. Proceedings of the International Society for Magnetic 
    Resonance in Medicine, 2010, 1670
    """
    input:
        fod=rules.mtnormalise.output.wm_fod,
        mask=rules.nii2mif.output.mask,
        convex_hull=rules.create_convex_hull.output.convex_hull,
        subcortical_seg=rules.get_num_nodes.input.seg,
    params:
        step=config["step"],
        sl=config["sl_count"],
    output:
        tck=bids_tractography_out(
            desc="iFOD2",
            suffix="tractography.tck",
        ),
    threads: 32
    resources:
        tmp_dir=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            **wildcards,
        ),
        tmp_tck=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            desc="iFOD2",
            suffix="tractography.tck",
            **wildcards,
        ),
        mem_mb=128000,
        time=60 * 24,
    log:
        bids_log(suffix="tckgen.log"),
    container:
        config["singularity"]["scattr"]
    conda:
        "../envs/mrtrix3.yaml"
    shell:
        """
        mkdir -p {resources.tmp_dir} 

        tckgen -nthreads {threads} -algorithm iFOD2 -step {params.step} \\
            -select {params.sl} -exclude {input.convex_hull} \\
            -include {input.subcortical_seg} -mask {input.mask} \\
            -seed_image {input.mask} {input.fod} {resources.tmp_tck} &> {log} 

        rsync -v {resources.tmp_tck} {output.tck} >> {log} 2>&1
        """


rule tcksift2:
    """
    Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The effects 
    of SIFT on the reproducibility and biological accuracy of the structural 
    connectome. NeuroImage, 2015, 104, 253-265
    """
    input:
        tck=rules.tckgen.output.tck,
        fod=rules.mtnormalise.output.wm_fod,
    output:
        weights=bids_tractography_out(
            desc="iFOD2",
            suffix="tckWeights.txt",
        ),
        mu=bids_tractography_out(desc="iFOD2", suffix="muCoefficient.txt"),
    threads: 8
    resources:
        mem_mb=32000,
        time=60 * 2,
    log:
        bids_log(suffix="tcksift2.log"),
    container:
        config["singularity"]["scattr"]
    conda:
        "../envs/mrtrix3.yaml"
    shell:
        """
        tcksift2 -nthreads {threads} -out_mu {output.mu} \\
            {input.tck} {input.fod} {output.weights} &> {log}
        """


checkpoint create_roi_mask:
    input:
        subcortical_seg=rules.tckgen.input.subcortical_seg,
        num_labels=rules.get_num_nodes.output.num_labels,
    params:
        base_dir=mrtrix_dir,
        subj_wildcards=inputs_dwi.subj_wildcards,
    output:
        out_dir=directory(bids_anat_out(datatype="roi_masks")),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    group:
        "tract_masks"
    container:
        config["singularity"]["scattr"]
    script:
        "../../scripts/mrtpipelines/create_roi_mask.py"


def aggregate_rois(wildcards):
    """Grab all created roi masks"""
    roi_masks = bids_anat_out(
        datatype="roi_masks",
        desc="{node}",
        suffix="mask.mif",
        **wildcards,
    )
    # Get node label wildcard
    node = glob_wildcards(roi_masks).node

    # Build node pairs
    node_pairs = np.triu_indices(len(node), k=1)

    return {
        "roi1": expand(
            bids_anat_out(
                datatype="roi_masks",
                desc="{node1}",
                suffix="mask.mif",
            ),
            node1=list(node_pairs[0] + 1),
            allow_missing=True,
        ),
        "roi2": expand(
            bids_anat_out(
                datatype="roi_masks",
                desc="{node2}",
                suffix="mask.mif",
            ),
            node2=list(node_pairs[1] + 1),
            allow_missing=True,
        ),
    }


checkpoint create_exclude_mask:
    input:
        unpack(aggregate_rois),
        rules.create_roi_mask.output.out_dir,
        subcortical_seg=rules.tckgen.input.subcortical_seg,
        num_labels=rules.get_num_nodes.output.num_labels,
    params:
        base_dir=mrtrix_dir,
        mask_dir=bids_anat_out(
            datatype="roi_masks",
        ),
        subj_wildcards=inputs_dwi.subj_wildcards,
    output:
        out_dir=directory(bids_anat_out(datatype="exclude_mask")),
    threads: 4
    resources:
        mem_mb=16000,
        time=60 * 3,
    group:
        "tract_masks"
    container:
        config["singularity"]["scattr"]
    script:
        "../../scripts/mrtpipelines/create_exclude_mask.py"


# TODO (v0.2): ADD OPTION TO OUTPUT TDI MAP


rule tck2connectome:
    """
    Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The 
    effects of SIFT on the reproducibility and biological accuracy of 
    the structural connectome. NeuroImage, 2015, 104, 253-265

    TODO (v0.2): ADD OPTION TO MULTIPLY BY MU COEFFICIENT
    """
    input:
        weights=rules.tcksift2.output.weights,
        tck=rules.tckgen.output.tck,
        subcortical_seg=rules.tckgen.input.subcortical_seg,
    params:
        radius=config["radial_search"],
    output:
        sl_assignment=bids_tractography_out(
            desc="subcortical",
            suffix="nodeAssignment.txt",
        ),
        node_weights=bids_tractography_out(
            desc="subcortical",
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
            desc="subcortical",
            suffix="nodeAssignment.txt",
            **wildcards,
        ),
        tmp_node_weights=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            desc="subcortical",
            suffix="nodeWeights.csv",
            **wildcards,
        ),
        mem_mb=128000,
        time=60 * 3,
    log:
        bids_log(suffix="tck2connectome.log"),
    group:
        "tractography_update"
    container:
        config["singularity"]["scattr"]
    conda:
        "../envs/mrtrix3.yaml"
    shell:
        """
        mkdir -p {resources.tmp_dir}

        tck2connectome -nthreads {threads} -zero_diagonal -stat_edge sum \\
        -assignment_radial_search {params.radius} \\
        -tck_weights_in {input.weights} \\
        -out_assignments {resources.tmp_sl_assignment} \\
        -symmetric {input.tck} {input.subcortical_seg} \\
        {resources.tmp_node_weights} &> {log}

        rsync {resources.tmp_sl_assignment} \\
        {output.sl_assignment} >> {log} 2>&1 

        rsync {resources.tmp_node_weights} {output.node_weights} >> {log} 2>&1
        """


checkpoint connectome2tck:
    input:
        node_weights=rules.tcksift2.output.weights,
        sl_assignment=rules.tck2connectome.output.sl_assignment,
        tck=rules.tckgen.output.tck,
        num_labels=rules.get_num_nodes.output.num_labels,
    output:
        output_dir=directory(
            str(
                Path(
                    bids_tractography_out(
                        datatype="unfiltered",
                    )
                ).parent
            )
        ),
    threads: 32
    resources:
        tmp_dir=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            datatype="unfiltered",
            **wildcards,
        ),
        edge_weight_prefix=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            datatype="unfiltered",
            desc="subcortical",
            suffix="tckWeights",
            **wildcards,
        ),
        edge_tck_prefix=lambda wildcards: bids_tractography_out(
            root=os.environ.get("SLURM_TMPDIR")
            if config.get("slurm_tmpdir")
            else "/tmp",
            datatype="unfiltered",
            desc="subcortical",
            suffix="from",
            **wildcards,
        ),
        mem_mb=128000,
        time=60 * 3,
    log:
        bids_log(suffix="connectome2tck.log"),
    group:
        "tractography_update"
    container:
        config["singularity"]["scattr"]
    conda:
        "../envs/mrtrix3.yaml"
    shell:
        """
        mkdir -p {resources.tmp_dir} {output.output_dir}

        num_labels=$(cat {input.num_labels})

        connectome2tck -nthreads {threads} \\
        -nodes `seq -s, 1 $((num_labels))` \\
        -exclusive -files per_edge -tck_weights_in {input.node_weights} \\
        -prefix_tck_weights_out {resources.edge_weight_prefix} \\
        {input.tck} {input.sl_assignment} \\
        {resources.edge_tck_prefix} &> {log}

        rsync {resources.edge_weight_prefix}*.csv \\
        {output.output_dir}/ >> {log} 2>&1
        rsync {resources.edge_tck_prefix}*.tck \\
        {output.output_dir}/ >> {log} 2>&1
        """

# Directories
zona_dir = str(Path(config["output_dir"]) / "zona_bb_subcortex")
labelmerge_dir = str(Path(config["output_dir"]) / "labelmerge")
log_dir = str(Path(config["output_dir"]) / ".logs" / "zona_bb_subcortex")

# Make directory if it doesn't exist
Path(zona_dir).mkdir(parents=True, exist_ok=True)


# BIDS partials
bids_anat = partial(
    bids,
    root=zona_dir,
    datatype="anat",
    **inputs.subj_wildcards,
)

bids_labelmerge = partial(
    bids,
    root=str(Path(labelmerge_dir) / "combined")
    if not config.get("skip_labelmerge")
    else config.get("labelmerge_base_dir") or zona_dir,
    **inputs.subj_wildcards,
)

bids_log = partial(
    bids,
    root=log_dir,
    **inputs.subj_wildcards,
)

""" References:
J.C. Lau, Y. Xiao, R.A.M. Haast, G. Gilmore, K. Uludağ, K.W. MacDougall, 
R.S. Menon, A.G. Parrent, T.M. Peters, A.R. Khan. Direct visualization and 
characterization of the human zona incerta and surrounding structures. 
Hum. Brain Mapp., 41 (2020), pp. 4500-4517, 10.1002/hbm.25137

Y. Xiao, J.C. Lau, T. Anderson, J. DeKraker, D.L. Collins, T. Peters, 
A.R. Khan. An accurate registration of the BigBrain dataset with the MNI PD25 
and ICBM152 atlases. Sci. Data, 6 (2019), p. 210, 10.1038/s41597-019-0217-0
"""


rule cp_zona_tsv:
    """Copy tsv to zona dir"""
    input:
        zona_tsv=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"]["tsv"])
        ),
    output:
        zona_tsv=f"{zona_dir}/desc-ZonaBB_dseg.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        time=10,
    group:
        "dseg_tsv"
    shell:
        "cp -v {input.zona_tsv} {output.zona_tsv}"


rule reg2native:
    """
    Create transforms from chosen template space to subject native space via 
    ANTsRegistrationSyNQuick
    """
    input:
        template=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"][config["Space"]]["dir"])
            / Path(config["zona_bb_subcortex"][config["Space"]]["T1w"])
        ),
        target=lambda wildcards: expand(
            inputs["T1w"].path,
            zip,
            **filter_list(inputs["T1w"].input_zip_lists, wildcards)
        )[0],
    params:
        out_dir=directory(str(Path(bids_anat()).parent)),
        out_prefix=bids_anat(
            desc=f"from{config['Space']}toNative_",
        ),
    output:
        t1w_nativespace=bids_anat(
            desc=f"from{config['Space']}toNative",
            suffix="Warped.nii.gz",
        ),
        warp=bids_anat(
            desc=f"from{config['Space']}toNative",
            suffix="1Warp.nii.gz",
        ),
        affine=bids_anat(
            desc=f"from{config['Space']}toNative",
            suffix="0GenericAffine.mat",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        bids_log(suffix="reg2native.log"),
    group:
        "subcortical_1"
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads}

        mkdir -p {params.out_dir}

        antsRegistrationSyNQuick.sh -n {threads} -d 3 \\ 
            -f {input.target} -m {input.template} \\
            -o {params.out_prefix} &> {log}
        """


rule warp2native:
    """Warp subcortical parcellations to subject native space"""
    input:
        dseg=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"][config["Space"]]["dir"])
            / Path(config["zona_bb_subcortex"][config["Space"]]["seg"])
        ),
        target=rules.reg2native.input.target,
        warp=rules.reg2native.output.warp,
        affine=rules.reg2native.output.affine,
    output:
        nii=bids_anat(
            space="T1w",
            desc="ZonaBB",
            suffix="dseg.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=30,
    log:
        bids_log(suffix="warp2native.log"),
    group:
        "subcortical_1"
    container:
        config["singularity"]["ants"]
    shell:
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} 

        antsApplyTransforms -v -d 3 -n MultiLabel \\
            -i {input.dseg} -r {input.target} \\
            -t {input.warp} -t {input.affine} \\
            -o {output.nii} &> {log}
        """


rule labelmerge:
    input:
        zona_seg=expand(
            rules.warp2native.output.nii
            if not config.get("labelmerge_base_dir")
            else [],
            zip,
            **subj_zip_list,
            allow_missing=True,
        ),
        fs_seg=expand(
            rules.fs_xfm_to_native.output.thal
            if not config.get("labelmerge_overlay_dir")
            else [],
            zip,
            **subj_zip_list,
            allow_missing=True,
        ),
        fs_tsv=rules.cp_fs_tsv.output.fs_tsv
        if not config.get("labelmerge_overlay_dir")
        else [],
        zona_tsv=rules.cp_zona_tsv.output.zona_tsv
        if not config.get("labelmerge_base_dir")
        else [],
    params:
        labelmerge_out_dir=directory(labelmerge_dir),
        labelmerge_base_dir=(
            config.get("labelmerge_base_dir")
            if config.get("labelmerge_base_dir")
            else zona_dir
        ),
        base_drops=f"--base_drops {config['labelmerge_base_drops']}",
        base_desc=f"--base_desc {config['labelmerge_base_desc']}",
        base_exceptions=(
            f"--base_exceptions {config.get('labelmerge_base_exceptions', '')}"
        ),
        overlay_dir=(
            f"--overlay_bids_dir {config.get('labelmerge_overlay_dir', rules.thalamic_segmentation.input.freesurfer_dir)}"
        ),
        overlay_desc=f"--overlay_desc {config['labelmerge_overlay_desc']}",
        overlay_drops=(
            f"--overlay_drops {config.get('labelmerge_overlay_drops', '')}"
        ),
        overlay_exceptions=(
            f"--overlay_exceptions {config.get('labelmerge_overlay_exceptions', '')}"
        ),
    output:
        seg=expand(
            bids_labelmerge(
                space="T1w",
                desc="combined",
                suffix="dseg.nii.gz",
            ),
            zip,
            **subj_zip_list,
            allow_missing=True,
        ),
        tsv=expand(
            bids_labelmerge(
                space="T1w",
                desc="combined",
                suffix="dseg.tsv",
            ),
            zip,
            **subj_zip_list,
            allow_missing=True,
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    group:
        "subcortical_group"
    container:
        config["singularity"]["labelmerge"]
    shell:
        """
        labelmerge {params.labelmerge_base_dir} {params.labelmerge_out_dir} \\
            participant \\ 
            {params.base_desc} {params.base_drops} {params.base_exceptions} \\
            {params.overlay_dir} {params.overlay_desc} \\ 
            {params.overlay_drops} {params.overlay_exceptions} \\
            --cores {threads} --force-output
        """


rule get_num_nodes:
    input:
        seg=bids_labelmerge(
            space="T1w",
            datatype="anat" if config.get("skip_labelmerge") else "",
            desc="combined"
            if not config.get("skip_labelmerge")
            else config.get("labelmerge_base_desc"),
            suffix="dseg.nii.gz",
        ),
    output:
        num_labels=temp(
            bids_labelmerge(
                space="T1w",
                datatype="anat" if config.get("skip_labelmerge") else "",
                desc="combined"
                if not config.get("skip_labelmerge")
                else config.get("labelmerge_base_desc"),
                suffix="numNodes.txt",
            )
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
    group:
        "subcortical_2"
    container:
        config["singularity"]["scattr"]
    script:
        "../scripts/zona_bb_subcortex/get_num_labels.py"


rule binarize:
    input:
        seg=rules.get_num_nodes.input.seg,
    output:
        mask=bids_labelmerge(
            space="T1w",
            desc="combined"
            if not config.get("skip_labelmerge")
            else config.get("labelmerge_base_desc"),
            suffix="mask.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
    log:
        bids_log(suffix="binarize.log"),
    group:
        "subcortical_2"
    container:
        config["singularity"]["scattr"]
    script:
        "../scripts/zona_bb_subcortex/create_seg_mask.py"


rule add_brainstem:
    input:
        mask=rules.binarize.output.mask,
        aparcaseg=rules.fs_xfm_to_native.output.aparcaseg
        if not config.get("skip_brainstem")
        else [],
    output:
        mask=bids_labelmerge(
            space="T1w",
            desc="labelmergeStem",
            suffix="mask.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
    log:
        bids_log(suffix="addBrainstem.log"),
    group:
        "subcortical_2"
    container:
        config["singularity"]["scattr"]
    script:
        "../scripts/zona_bb_subcortex/add_brainstem.py"


rule create_convex_hull:
    input:
        bin_seg=rules.binarize.output.mask
        if config["skip_brainstem"]
        else rules.add_brainstem.output.mask,
    output:
        convex_hull=bids_labelmerge(
            space="T1w",
            desc="ConvexHull",
            suffix="mask.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        bids_log(suffix="createConvexHull.log"),
    group:
        "subcortical_2"
    container:
        config["singularity"]["scattr"]
    script:
        "../scripts/zona_bb_subcortex/convexHull_roi.py"

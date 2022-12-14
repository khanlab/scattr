# Directories
zona_dir = str(Path(config["output_dir"]) / "zona_bb_subcortex")
labelmerge_dir = str(Path(config["output_dir"]) / "labelmerge")

# Make directory if it doesn't exist
Path(zona_dir).mkdir(parents=True, exist_ok=True)


# BIDS partials
bids_anat = partial(
    bids,
    root=zona_dir,
    datatype="anat",
    **config["subj_wildcards"],
)

bids_labelmerge = partial(
    bids,
    root=str(Path(labelmerge_dir) / "combined"),
    **config["subj_wildcards"],
)

# References:
# J.C. Lau, Y. Xiao, R.A.M. Haast, G. Gilmore, K. Uludağ, K.W. MacDougall, R.S. Menon, A.G. Parrent, T.M. Peters, A.R. Khan. Direct visualization and characterization of the human zona incerta and surrounding structures. Hum. Brain Mapp., 41 (2020), pp. 4500-4517, 10.1002/hbm.25137

# Y. Xiao, J.C. Lau, T. Anderson, J. DeKraker, D.L. Collins, T. Peters, A.R. Khan. An accurate registration of the BigBrain dataset with the MNI PD25 and ICBM152 atlases. Sci. Data, 6 (2019), p. 210, 10.1038/s41597-019-0217-0


rule cp_zona_tsv:
    """Copy tsv to zona dir"""
    input:
        zona_tsv=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"]["tsv"])
        ),
    output:
        zona_tsv=f"{zona_dir}/desc-ZonaBBSubcor_dseg.tsv",
    threads: 4
    resources:
        mem_mb=8000,
        time=10,
    group: "dseg_tsv"
    shell:
        "cp -v {input.zona_tsv} {output.zona_tsv}"


rule reg2native:
    """Create transfroms from chosen template space to subject native space via ANTsRegistrationSyNQuick"""
    input:
        template=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"][config["Space"]]["dir"])
            / Path(config["zona_bb_subcortex"][config["Space"]]["T1w"])
        ),
        target=config["input_path"]["T1w"], 
    params:
        out_dir=directory(bids_anat()),
        out_prefix=bids_anat(
            desc=f"from{config['Space']}toNative_",
        ),
    output:
        warp = bids_anat(
            desc=f"from{config['Space']}toNative",
            suffix="1Warp.nii.gz",
        ),
        affine = bids_anat(
            desc=f"from{config['Space']}toNative",
            suffix="0GenericAffine.mat",
        ),
    threads: 8
    resources:
        mem_mb=32000,
        time=60,
    log:
        f"{config['output_dir']}/logs/zona_bb_subcortex/sub-{{subject}}/reg2native.log",
    group: "subcortical"
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "mkdir -p {params.out_dir} && "
        "export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} && "
        "antsRegistrationSyNQuick.sh -n {threads} -d 3 " 
        "-f {input.target} -m {input.template} "
        "-o {params.out_prefix} &> {log}"

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
        )
    threads: 8
    resources:
        mem_mb=32000,
        time=10,
    log:
        f"{config['output_dir']}/logs/zona_bb_subcortex/sub-{{subject}}/warp2native.log",
    group: "subcortical"
    container:
        config["singularity"]["ants"]
    shell:
        "export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} && "
        "antsApplyTransforms -v -d 3 -n MultiLabel "
        "-i {input.dseg} -r {input.target} "
        "-t {input.warp} -t {input.affine} "
        "-o {output.nii} &> {log}" 


# TODO: Add back later on
# rule xfm_zona_rois:
#     input:
#         mask=str(
#             Path(config["zona_bb_subcortex"][config["Space"]]["dir"])
#             / f'sub-SNSX32Nlin2020Asym_space-{config["Space"]}_hemi-{{hemi}}_desc-{{struct}}_mask.nii.gz'
#         ),
#         ref=rules.xfm2native.input.ref,
#         xfm=rules.xfm2native.output.xfm,
#     output:
#         mask=bids(
#             root=zona_dir,
#             datatype="anat",
#             space="T1w",
#             hemi="{{hemi,(L|R)}}",
#             desc="{{struct,(fct|ft|fl|hfields)}}",
#             suffix="mask.nii.gz",
#         ),
#     container:
#         config["singularity"]["neuroglia-core"]
#     shell:
#         "applywarp --rel --interp=nn -i {input.mask} -r {input.ref} -w {input.xfm} -o {output.mask}"


rule rm_bb_thal:
    """Removes thalamus from existing parcellation

    1. Grab labels following thalamus ROI
    2. Remove existing thalamus from segmentation
    3. Add labels following thalamus ROI back to segmentation
    """
    input:
        seg=rules.warp2native.output.nii,
    output:
        seg=bids_anat(
            space="T1w",
            desc="ZonaBBSubcor",
            suffix="dseg.nii.gz",
        ),
        non_thal=temp(
            bids_anat(
                space="T1w",
                desc="NonThal",
                suffix="dseg.nii.gz",
            )
        ),
        rm_seg=temp(
            bids_anat(
                space="T1w",
                desc="ThalPost",
                suffix="dseg.nii.gz",
            )
        ),
    threads: 8
    resources:
        mem_mb=32000,
        time=10,
    log:
        f"{config['output_dir']}/logs/zona_bb_subcortex/sub-{{subject}}/rm_bb_thal.log",
    group: "subcortical"
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.seg} -sub 2 {output.non_thal} &> {log} && "
        "fslmaths {output.non_thal} -thr 15 {output.non_thal} >> {log} 2>&1 && "
        "fslmaths {input.seg} -thr 15 {output.rm_seg} >> {log} 2>&1 && "
        "fslmaths {input.seg} -sub {output.rm_seg} {output.seg} >> {log} 2>&1 && "
        "fslmaths {output.seg} -add {output.non_thal} {output.seg} >> {log} 2>&1"


rule labelmerge:
    input:
        zona_seg=expand(
            bids_anat(
                subject="{subject}",
                space="T1w",
                desc="ZonaBBSubcor",
                suffix="dseg.nii.gz",
            ),
            subject=config["input_lists"]["T1w"]["subject"],
        ),
        fs_seg=expand(
            rules.fs_xfm_to_native.output.thal, 
            subject=config["input_lists"]["T1w"]["subject"],
        )
    params:
        zona_dir=zona_dir,
        fs_dir=rules.thalamic_segmentation.input.freesurfer_dir,
        zona_desc="ZonaBBSubcor",
        fs_desc="FreesurferThal",
        labelmerge_dir=directory(labelmerge_dir),
        labelmerge_container=config["singularity"]["labelmerge"],
    output:
        seg=expand(
            bids_labelmerge(
                subject="{subject}",
                space="T1w",
                desc="combined",
                suffix="dseg.nii.gz",
            ),
            subject=config["input_lists"]["T1w"]["subject"],
        ),
        tsv=expand(
            bids_labelmerge(
                subject="{subject}",
                space="T1w",
                desc="combined",
                suffix="dseg.tsv",
            ),
            subject=config["input_lists"]["T1w"]["subject"],
        ),
    threads: 8
    resources:
        mem_mb=32000,
        time=10,
    group: "subcortical"
    shell:
        "singularity run {params.labelmerge_container} {params.zona_dir} {params.labelmerge_dir} participant --base_desc {params.zona_desc} --overlay_bids_dir {params.fs_dir} --overlay_desc {params.fs_desc} -c all --force-output"


rule get_num_nodes:
    input:
        seg=rules.labelmerge.output.seg,
    output:
        num_labels=temp(
            bids_labelmerge(
                space="T1w",
                desc="combined",
                suffix="numNodes.txt",
            )
        )
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
    group: "subcortical"
    container:
        config["singularity"]["dbsc"]
    script:
        "../scripts/zona_bb_subcortex/get_num_labels.py"

rule binarize:
    input:
        seg=bids_labelmerge(
            space="T1w",
            desc="combined",
            suffix="dseg.nii.gz",
        ),
    output:
        mask=bids_labelmerge(
            space="T1w",
            desc="combined",
            suffix="mask.nii.gz",
        ),
    threads: 8
    resources:
        mem_mb=16000,
        time=10,
    log:
        f"{config['output_dir']}/logs/labelmerge/sub-{{subject}}/binarize.log",
    group: "subcortical"
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.seg} -bin {output.mask} &> {log}"


rule add_brainstem:
    input:
        mask=rules.binarize.output.mask,
        aparcaseg=rules.fs_xfm_to_native.output.aparcaseg,
    output:
        mask=bids_anat(
            root=labelmerge_dir,
            space="T1w",
            desc="labelmergeStem",
            suffix="mask.nii.gz",
        ),
    threads: 8
    resources:
        mem_mb=32000,
        time=10,
    log:
        f"{config['output_dir']}/logs/labelmerge/sub-{{subject}}/add_brainstem.log",
    group: "subcortical"
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.aparcaseg} -thr 16 -uthr 16 -bin -max {input.mask} {output.mask} &> {log}"


rule create_convex_hull:
    input:
        bin_seg=rules.add_brainstem.output.mask,
    output:
        convex_hull=bids_anat(
            root=labelmerge_dir,
            space="T1w",
            desc="ConvexHull",
            suffix="mask.nii.gz",
        ),
    threads: 8
    resources:
        mem_mb=32000,
        time=30,
    log:
        f"{config['output_dir']}/logs/labelmerge/sub-{{subject}}/create_convex_hull.log",
    group: "subcortical"
    container:
        config["singularity"]["dbsc"]
    script:
        "../scripts/zona_bb_subcortex/convexHull_roi.py"

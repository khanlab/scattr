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
    root=labelmerge_dir,
    datatype="anat",
    **config["subj_wildcards"],
)

# References:
# J.C. Lau, Y. Xiao, R.A.M. Haast, G. Gilmore, K. Uludağ, K.W. MacDougall, R.S. Menon, A.G. Parrent, T.M. Peters, A.R. Khan. Direct visualization and characterization of the human zona incerta and surrounding structures. Hum. Brain Mapp., 41 (2020), pp. 4500-4517, 10.1002/hbm.25137

# Y. Xiao, J.C. Lau, T. Anderson, J. DeKraker, D.L. Collins, T. Peters, A.R. Khan. An accurate registration of the BigBrain dataset with the MNI PD25 and ICBM152 atlases. Sci. Data, 6 (2019), p. 210, 10.1038/s41597-019-0217-0


rule cp_zona_tsv:
    """Copy tsv to zona dir"""
    input:
        zona_tsv=str(
            Path(workflow.basedir).parent / Path(config["zona_bb_subcortex"]["tsv"])
        ),
    output:
        zona_tsv=f"{zona_dir}/desc-ZonaBBSubcor_dseg.tsv",
    shell:
        "cp -v {input.zona_tsv} {output.zona_tsv}"


rule xfm2native:
    """Transform from subcortical parcellations from chosen template space to subject native space"""
    input:
        seg=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"][config["Space"]]["dir"])
            / Path(config["zona_bb_subcortex"][config["Space"]]["seg"])
        ),
        seg_anat=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"][config["Space"]]["dir"])
            / Path(config["zona_bb_subcortex"][config["Space"]]["T1w"])
        ),
        ref=config["input_path"]["T1w"],
    output:
        xfm=bids_anat(
            desc=f"from{config['Space']}toNative",
            suffix="xfm.mat",
        ),
        nii=bids_anat(
            space="T1w",
            desc="ZonaBB",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "flirt -in {input.seg_anat} -r {input.ref} -omat {output.xfm} && "
        "applywarp --rel --interp=nn -i {input.seg} -r {input.ref} -w {output.xfm} -o {output.nii}"


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
        seg=rules.xfm2native.output.nii,
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
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.seg} -sub 2 {output.non_thal} && "
        "fslmaths {output.non_thal} -thr 15 {output.non_thal} && "
        "fslmaths {input.seg} -thr 15 {output.rm_seg} && "
        "fslmaths {input.seg} -sub {output.rm_seg} {output.seg} && "
        "fslmaths {output.seg} -add {output.non_thal} {output.seg}"


rule labelmerge:
    input:
        seg=expand(
            bids_anat(
                subject="{subject}",
                space="T1w",
                desc="ZonaBBSubcor",
                suffix="dseg.nii.gz",
            ),
            subject=config["input_lists"]["dwi"]["subject"],
        ),
    params:
        zona_dir=zona_dir,
        fs_dir=str(Path(config["output_dir"]) / "freesurfer"),
        zona_desc="ZonaBBSubcor",
        fs_desc="FreesurferThal",
        labelmerge_dir=directory(labelmerge_dir),
    output:
        seg=expand(
            bids_labelmerge(
                subject="{subject}",
                space="T1w",
                desc="combined",
                suffix="dseg.nii.gz",
            ),
            subject=config["input_lists"]["dwi"]["subject"],
        ),
        tsv=expand(
            bids_labelmerge(
                subject="{subject}",
                space="T1w",
                desc="combined",
                suffix="dseg.tsv",
            ),
            subject=config["input_lists"]["dwi"]["subject"],
        ),
    # ADD CONTAINER
    shell:
        # TO BE UPDATED WITH APPROPRIATE COMMAND
        "run.py {params.zona_dir} {params.labelmerge_dir} --overlay_bids_dir {params.fs_dir} --overlay_desc {params.fs_desc} --base_desc {params.zona_desc}"


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
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.seg} -bin {output.mask}"


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
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.aparcaseg} -thr 16 -uthr 16 -bin -max {input.mask} {output.mask}"


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
    script:
        "../resources/zona_bb_subcortex/convexHull_roi.py"

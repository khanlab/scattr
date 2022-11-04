from os.path import join


# Rules to include
include: "freesurfer.smk"


# References:
# J.C. Lau, Y. Xiao, R.A.M. Haast, G. Gilmore, K. Uludağ, K.W. MacDougall, R.S. Menon, A.G. Parrent, T.M. Peters, A.R. Khan. Direct visualization and characterization of the human zona incerta and surrounding structures. Hum. Brain Mapp., 41 (2020), pp. 4500-4517, 10.1002/hbm.25137

# Y. Xiao, J.C. Lau, T. Anderson, J. DeKraker, D.L. Collins, T. Peters, A.R. Khan. An accurate registration of the BigBrain dataset with the MNI PD25 and ICBM152 atlases. Sci. Data, 6 (2019), p. 210, 10.1038/s41597-019-0217-0

# Directories
zona_dir = join(config["output_dir"], "zona_bb_subcortex")


rule xfm2native:
    """Transform from Zona template space to subject native space"""
    input:
        seg=config["zona_bb_subcortex"][config["Space"]]["seg"],
        seg_anat=config["zona_bb_subcortex"][config["Space"]]["T1w"],
        ref=bids(
            root=config["bids_dir"],
            datatype="anat",
            suffix="T1w.nii.gz",
            **config["subj_wildcards"],
        ),
    output:
        xfm=bids(
            root=zona_dir,
            datatype="anat",
            desc=f"from{config['Space']}toNative",
            suffix="xfm.mat",
            **config["subj_wildcards"],
        ),
        nii=bids(
            root=zona_dir,
            datatype="anat",
            space="T1w",
            desc="ZonaBB",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "flirt -in {input.seg_anat} -r {input.ref} -omat {output.xfm} &&"
        "applywarp --rel --interp=nn -i {input.seg} -r {input.ref} -w {output.xfm} -o {output.nii}"


rule binarize:
    input:
        nii=rules.xfm2native.output.nii,
    output:
        bin=bids(
            root=zona_dir,
            datatype="anat",
            space="T1w",
            desc="ZonaBB",
            suffix="mask.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.nii} -bin {output.bin}"


rule add_brainstem:
    input:
        binarize=rules.binarize.output.bin,
        aparcaseg=rules.fs_xfm2native.output.aparcaseg,
    output:
        binarize=bids(
            root=zona_dir,
            datatype="anat",
            space="T1w",
            desc="ZonaBBStem",
            suffix="mask.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.aparcaseg} -thr 16 -uthr 16 -bin -max {input.binarize} {output.binarize}"


rule xfm_zona_rois:
    input:
        mask=expand(
            join(
                config["zona_bb_subcortex"][config["Space"]]["dir"],
                f'sub-SNSX32Nlin2020Asym_space-{config["Space"]}_hemi-{{hemi}}_desc-{{struct}}_mask.nii.gz',
            ),
            hemi=["L", "R"],
            struct=["fct", "ft", "fl", "hfields"],
        ),
        ref=rules.xfm2native.input.ref,
        xfm=rules.xfm2native.input.xfm,
    output:
        mask=expand(
            bids(
                root=zona_dir,
                datatype="anat",
                space="T1w",
                hemi="{hemi}",
                desc="{struct}",
                suffix="mask.nii.gz",
            ),
            hemi=["L", "R"],
            struct=["fct", "ft", "fl", "hfields"],
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "applywarp --rel --interp=nn -i {input.mask} -r {input.ref} -w {input.xfm} -o {output.mask}"


# Everything in the block below should be replaced by labelmerge
######################################################################
rule rm_bb_thal:
    """Removes existing thalamus"""
    input:
        seg=rules.xfm2native.output.nii,
    output:
        seg=bids(
            root=join(config["output_dir"], "zona_bb_subcortex"),
            datatype="anat",
            space="T1w",
            desc="ZonaBBNoThal",
            suffix="seg.nii.gz",
            **config["subj_wildcards"],
        ),
        non_thal=temp(
            bids(
                root=join(config["output_dir"], "zona_bb_subcortex"),
                datatype="anat",
                space="T1w",
                desc="NonThal",
                suffix="seg.nii.gz",
                **config["subj_wildcards"],
            )
        ),
        rm_seg=temp(
            bids(
                root=join(config["output_dir"], "zona_bb_subcortex"),
                datatype="anat",
                space="T1w",
                desc="rm",
                suffix="seg.nii.gz",
                **config["subj_wildcards"],
            )
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.seg} -sub 2 {output.non_thal} "
        "fslmaths {output.non_thal} -thr 15 {output.non_thal} "
        "fslmaths {input.seg} -thr 15 {output.rm_seg} "
        "fslmaths {input.seg} -sub {output.rm_seg} {output.seg}"
        "fslmaths {output.seg} -add {output.non_thal} {output.seg}"


rule add_new_thal:
    input:
        aparcaseg=rules.add_brainstem.input.aparcaseg,
        labels=config["freesurfer"]["labels"],
        thal=rules.fs_xfm2native.output.thal,
        seg=rules.rm_bb_thal.output.rm_seg,
    output:
        seg=bids(
            root=join(config["output_dir"], "zona_bb_subcortex"),
            datatype="anat",
            space="T1w",
            desc="ZonaBBFSThal",
            suffix="seg.nii.gz",
            **config["subj_wildcards"],
        ),
        label=temp(
            bids(
                root=join(config["output_dir"], "zona_bb_subcortex"),
                datatype="anat",
                space="T1w",
                suffix="label.nii.gz",
                **config["subj_wildcards"],
            )
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "labelIdx=23 "
        "cp {input.seg} {output.seg} "
        'while IFS=, read -r label hemi nuclei || [ -n "$label" ]; do '
        '  if [ "$label" = "# Label" ]; then '
        "    continue "
        "  fi "
        "  fslmaths {input.thal} -thr $label -uthr $label -bin {output.label} "
        "  fslmaths {output.label} -mul $labelIdx {output.label} "
        "  fslmaths {output.seg} -max {output.label} {output.seg} "
        "  labelIdx=$((labelIdx+1)) "
        " done < {input.labels}"


rule bin_new_seg:
    input:
        seg=rules.add_new_thal.output.seg,
    output:
        seg=bids(
            root=join(config["output_dir"], "zona_bb_subcortex"),
            datatype="anat",
            space="T1w",
            desc="ZonaBBFSThal",
            suffix="seg_bin.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.seg} -bin {output.seg}"


rule add_brainstem_new_seg:
    input:
        seg=rules.bin_new_seg.output.seg,
        aparcaseg=rules.add_brainstem.input.aparcaseg,
    output:
        seg=bids(
            root=join(config["output_dir"], "zona_bb_subcortex"),
            datatype="anat",
            space="T1w",
            desc="ZonaBBFSThalStem",
            suffix="seg_bin.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.aparcaseg} -thr 16 -uthr 16 -bin -max {input.seg} {output.seg}"


###################BLOCK TO REPLACE WITH LABELMERGE############################


rule create_convex_hull:
    input:
        bin_seg=rules.add_brainstem_new_seg.output.seg,
    output:
        convex_hull=bids(
            root=zona_dir,
            datatype="anat",
            space="T1w",
            desc="ConvexHull",
            suffix="mask.nii.gz",
            **config["subj_wildcards"],
        ),
    script:
        "../resources/zona_bb_subcortex/convexHull_roi.py"

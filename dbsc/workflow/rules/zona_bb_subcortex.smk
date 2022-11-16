from pathlib import Path
from functools import partial


# Directories
zona_dir = str(Path(config["output_dir"]) / "zona_bb_subcortex")

# Make directory if it doesn't exist
Path(zona_dir).mkdir(parents=True, exist_ok=True)


# BIDS partials
bids_anat = partial(
    bids, root=zona_dir, datatype="anat", **config["subj_wildcards"]
)

# References:
# J.C. Lau, Y. Xiao, R.A.M. Haast, G. Gilmore, K. Uludağ, K.W. MacDougall, R.S. Menon, A.G. Parrent, T.M. Peters, A.R. Khan. Direct visualization and characterization of the human zona incerta and surrounding structures. Hum. Brain Mapp., 41 (2020), pp. 4500-4517, 10.1002/hbm.25137

# Y. Xiao, J.C. Lau, T. Anderson, J. DeKraker, D.L. Collins, T. Peters, A.R. Khan. An accurate registration of the BigBrain dataset with the MNI PD25 and ICBM152 atlases. Sci. Data, 6 (2019), p. 210, 10.1038/s41597-019-0217-0


rule xfm2native:
    """Transform from Zona template space to subject native space"""
    input:
        seg=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"][config["Space"]]["seg"])
        ),
        seg_anat=str(
            Path(workflow.basedir).parent
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


rule binarize:
    input:
        nii=rules.xfm2native.output.nii,
    output:
        mask=bids_anat(
            space="T1w",
            desc="ZonaBB",
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.nii} -bin {output.mask}"


rule add_brainstem:
    input:
        binarize=rules.binarize.output.mask,
        aparcaseg=rules.fs_xfm_to_native.output.aparcaseg,
    output:
        binarize=bids_anat(
            space="T1w",
            desc="ZonaBBStem",
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "fslmaths {input.aparcaseg} -thr 16 -uthr 16 -bin -max {input.binarize} {output.binarize}"


rule xfm_zona_rois:
    input:
        mask=str(
            Path(config["zona_bb_subcortex"][config["Space"]]["dir"])
            / f'sub-SNSX32Nlin2020Asym_space-{config["Space"]}_hemi-{{hemi}}_desc-{{struct}}_mask.nii.gz'
        ),
        ref=rules.xfm2native.input.ref,
        xfm=rules.xfm2native.output.xfm,
    output:
        mask=bids(
            root=zona_dir,
            datatype="anat",
            space="T1w",
            hemi="{{hemi,(L|R)}}",
            desc="{{struct,(fct|ft|fl|hfields)}}",
            suffix="mask.nii.gz",
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
        seg=rules.fs_xfm_to_native.output.thal,
    output:
        seg=bids(
            root=zona_dir,
            datatype="anat",
            space="T1w",
            desc="ZonaBBNoThal",
            suffix="seg.nii.gz",
            **config["subj_wildcards"],
        ),
        non_thal=temp(
            bids(
                root=zona_dir,
                datatype="anat",
                space="T1w",
                desc="NonThal",
                suffix="seg.nii.gz",
                **config["subj_wildcards"],
            )
        ),
        rm_seg=temp(
            bids(
                root=zona_dir,
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
        labels=str(
            Path(workflow.basedir).parent
            / Path(config["freesurfer"]["labels"])
        ),
        thal=rules.fs_xfm_to_native.output.thal,
        seg=rules.rm_bb_thal.output.rm_seg,
    output:
        seg=bids(
            root=zona_dir,
            datatype="anat",
            space="T1w",
            desc="ZonaBBFSThal",
            suffix="seg.nii.gz",
            **config["subj_wildcards"],
        ),
        label=temp(
            bids(
                root=zona_dir,
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
            root=zona_dir,
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
            root=zona_dir,
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
        convex_hull=bids_anat(
            space="T1w",
            desc="ConvexHull",
            suffix="mask.nii.gz",
        ),
    script:
        "../resources/zona_bb_subcortex/convexHull_roi.py"

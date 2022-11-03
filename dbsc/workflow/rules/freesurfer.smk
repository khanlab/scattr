from os.path import join

# Directories
freesurfer_dir = join(config["out_dir"], "freesurfer")


rule thalamic_segmentation:
    """Perform thalamus segmentation

    NOTE: Outputname is defined by Freesurfer
    """
    input:
        freesurfer_dir=freesurfer_dir,
    params:
        fs_license=config["fs_license"],
    output:
        thal_seg=bids(
            root=freesurfer_dir,
            datatype="anat",
            suffix="T1.mgz",
            space="Freesurfer",
            **config["subj_wildcards"],
        ),
    threads: workflow.cores
    container:
        config["singularity"]["freesurfer"]
    shell:
        "FS_LICENSE={params.fs_license} "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={params.threads} "
        "SUBJECTS_DIR={input.freesurfer_dir} "
        "segmentThalamicNuclei.sh {{subject}}"


fs_out = bids(
    root=freesurfer_dir,
    datatype="anat",
    suffix="{suffix}",
    space="{space}",
    **config["subj_wildcards"],
)


rule mgz2nii:
    """Convert from .mgz to .nii.gz"""
    input:
        thal=rules.thalamic_segmentation.output.thal_seg,
        aparcaseg=join(freesurfer_dir, "{subject}/mri/aparc+aseg.mgz"),
        lRibbon=join(freesurfer_dir, "{subject}/mri/lh.ribbon.mgz"),
        rRibbon=join(freesurfer_dir, "{subject}/mri/rh.ribbon.mgz"),
    params:
        fs_license=config["fs_license"],
    output:
        thal=expand(
            fs_out, suffix="thalamus.nii.gz", space="Freesurfer", allow_missing=True
        ),
        aparcaseg=expand(
            fs_out, suffix="aparcaseg.nii.gz", space="Freesurfer", allow_missing=True
        ),
        ribbon_mgz=expand(
            temp(fs_out), suffix="ribbon.mgz", space="Freesurfer", allow_missing=True
        ),
        ribbon=expand(
            fs_out, suffix="ribbon.nii.gz", space="Freesurfer", allow_missing=True
        ),
    threads: workflow.cores
    container:
        config["singularity"]["freesurfer"]
    shell:
        "FS_LICENSE={params.fs_license} "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={params.threads} "
        "mri_convert {input.thal} {output.thal} "
        "mri_convert {input.aparcaseg} {output.aparcaseg} "
        "mergeseg --src {input.lRibbon} --merge {input.rRibbon} --o {output.ribbon_mgz}"
        "mri_convert {output.ribbon_mgz} {output.ribbon}"


rule fs_xfm_to_native:
    """Transform from Freesurfer space to T1w space"""
    input:
        thal=rules.mgz2nii.output.thal,
        aparcaseg=rules.mgz2nii.output.aparcaseg,
        ribbon=rules.mgz2nii.output.ribbon,
        ref=bids(
            root=config["bids_dir"],
            datatype="anat",
            suffix="T1w.nii.gz",
            **config["subj_wildcards"],
        ),
    output:
        thal=expand(fs_out, suffix="thalamus.nii.gz", space="T1w", allow_missing=True),
        aparcaseg=expand(
            fs_out, suffix="aparcaseg.nii.gz", space="T1w", allow_missing=True
        ),
        ribbon=expand(fs_out, suffix="ribbon.nii.gz", space="T1w", allow_missing=True),
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "antsApplyTransforms -d 3 -n MultiLabel -i {input.thal} -r {input.ref} -o {output.thal} "
        "antsApplyTransforms -d 3 -n MultiLabel -i {input.aparcaseg} -r {input.ref} -o {output.aparcaseg} "
        "antsApplyTransforms -d 3 -n MultiLabel -i {input.ribbon} -r {input.ref} -o {output.ribbon} "

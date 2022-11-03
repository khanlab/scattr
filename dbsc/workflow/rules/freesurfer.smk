from os.path import join

# Directories
freesurfer_dir = join(config["out_dir"], "freesurfer")

# Freesurfer references (with additional in rules as necessary)
# B. Fischl, A. van der Kouwe, C. Destrieux, E. Halgren, F. Ségonne, D.H. Salat, E. Busa, L.J. Seidman, J. Goldstein, D. Kennedy, V. Caviness, N. Makris, B. Rosen, A.M. Dale. Automatically parcellating the human cerebral cortex. Cereb. Cortex, 14 (2004), pp. 11-22, 10.1093/cercor/bhg087


rule thalamic_segmentation:
    """Perform thalamus segmentation

    Reference: J.E. Iglesias, R. Insausti, G. Lerma-Usabiaga, M. Bocchetta, K. Van Leemput, D.N. Greve, A. van der Kouwe, B. Fischl, C. Caballero-Gaudes, P.M. Paz-Alonso. A probabilistic atlas of the human thalamic nuclei combining ex vivo MRI and histology. NeuroImage, 183 (2018), pp. 314-326, 10.1016/j.neuroimage.2018.08.012

    NOTE: Output name is defined by Freesurfer
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
            fs_out, suffix="ribbon.mgz", space="Freesurfer", allow_missing=True
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

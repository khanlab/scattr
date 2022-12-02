# Directories
freesurfer_dir = str(Path(config["output_dir"]) / "freesurfer")
if config.get("freesurfer_dir"):
    freesurfer_dir = config.get("freesurfer_dir")

# BIDS partials
bids_fs_out = partial(
    bids,
    root=freesurfer_dir,
    datatype="anat",
    **config["subj_wildcards"],
)

# Freesurfer references (with additional in rules as necessary)
# B. Fischl, A. van der Kouwe, C. Destrieux, E. Halgren, F. Ségonne, D.H. Salat, E. Busa, L.J. Seidman, J. Goldstein, D. Kennedy, V. Caviness, N. Makris, B. Rosen, A.M. Dale. Automatically parcellating the human cerebral cortex. Cereb. Cortex, 14 (2004), pp. 11-22, 10.1093/cercor/bhg087


rule cp_fs_tsv:
    """Copy tsv to freesurfer dir
    """
    input:
        fs_tsv=str(
            Path(workflow.basedir).parent / Path(config["freesurfer"]["tsv"])
        ),
    output:
        fs_tsv=expand(
            "{freesurfer_dir}/desc-FreesurferThal_dseg.tsv",
            freesurfer_dir=freesurfer_dir,
        ),
    threads: 4
    resources:
        mem_mb=8000,
        time=10,
    shell:
        "cp -v {input.fs_tsv} {output.fs_tsv}"


rule thalamic_segmentation:
    """Perform thalamus segmentation

    Reference: J.E. Iglesias, R. Insausti, G. Lerma-Usabiaga, M. Bocchetta, K. Van Leemput, D.N. Greve, A. van der Kouwe, B. Fischl, C. Caballero-Gaudes, P.M. Paz-Alonso. A probabilistic atlas of the human thalamic nuclei combining ex vivo MRI and histology. NeuroImage, 183 (2018), pp. 314-326, 10.1016/j.neuroimage.2018.08.012

    NOTE: Output file is defined by Freesurfer script
    """
    input:
        freesurfer_dir=freesurfer_dir,
    params:
        fs_license=config["fs_license"],
    output:
        thal_seg=str(
            Path(freesurfer_dir)
            / "sub-{subject}/mri/ThalamicNuclei.v12.T1.mgz"
        ),
    threads: 8
    resources:
        mem_mb=16000,
        time=30,
    log:
        f"{config['output_dir']}/logs/freesurfer/sub-{{subject}}/thalamic_segmentation.log",
    group: "freesurfer"
    container:
        config["singularity"]["freesurfer"]
    shell:
        "FS_LICENSE={params.fs_license} && "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} && "
        "segmentThalamicNuclei.sh sub-{wildcards.subject} {input.freesurfer_dir} &> {log}"


rule mgz2nii:
    """
    Convert from .mgz to .nii.gz

    NOTE: During conversion, files are renamed to BIDS-esque formatting
    """
    input:
        thal=rules.thalamic_segmentation.output.thal_seg,
        aparcaseg=str(
            Path(freesurfer_dir) / "sub-{subject}/mri/aparc+aseg.mgz"
        ),
        lRibbon=str(Path(freesurfer_dir) / "sub-{subject}/mri/lh.ribbon.mgz"),
        rRibbon=str(Path(freesurfer_dir) / "sub-{subject}/mri/rh.ribbon.mgz"),
    params:
        fs_license=config["fs_license"],
    output:
        thal=bids_fs_out(space="Freesurfer", suffix="thalamus.nii.gz"),
        aparcaseg=bids_fs_out(space="Freesurfer", suffix="aparcaseg.nii.gz"),
        ribbon_mgz=bids_fs_out(space="Freesurfer", suffix="ribbon.mgz"),
        ribbon=bids_fs_out(space="Freesurfer", suffix="ribbon.nii.gz"),
    threads: 8
    resources:
        mem_mb=16000,
        time=10,
    log:
        f"{config['output_dir']}/logs/freesurfer/{{subject}}/mgz2nii.log",
    group: "freesurfer"
    container:
        config["singularity"]["freesurfer"]
    shell:
        "FS_LICENSE={params.fs_license} && "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} && "
        "mri_convert {input.thal} {output.thal} &> {log} &&  "
        "mri_convert {input.aparcaseg} {output.aparcaseg} >> {log} 2>&1 && "
        "mergeseg --src {input.lRibbon} --merge {input.rRibbon} --o {output.ribbon_mgz} >> {log} 2>&1 && "
        "mri_convert {output.ribbon_mgz} {output.ribbon} >> {log} 2>&1"


rule fs_xfm_to_native:
    """Transform from Freesurfer space to T1w space"""
    input:
        thal=rules.mgz2nii.output.thal,
        aparcaseg=rules.mgz2nii.output.aparcaseg,
        ribbon=rules.mgz2nii.output.ribbon,
        ref=config["input_path"]["T1w"],
    output:
        thal=bids_fs_out(space="T1w", suffix="thalamus.nii.gz"),
        aparcaseg=bids_fs_out(space="T1w", suffix="aparcaseg.nii.gz"),
        ribbon=bids_fs_out(space="T1w", suffix="ribbon.nii.gz"),
    threads: 8
    resources:
        mem_mb=16000,
        time=30,
    log:
        f"{config['output_dir']}/logs/freesurfer/{{subject}}/fs_xfm_to_native.log",
    group: "freesurfer"
    container:
        config["singularity"]["neuroglia-core"]
    shell:
        "antsApplyTransforms -d 3 -n MultiLabel -i {input.thal} -r {input.ref} -o {output.thal} &> {log} && "
        "antsApplyTransforms -d 3 -n MultiLabel -i {input.aparcaseg} -r {input.ref} -o {output.aparcaseg} >> {log} 2>&1 && "
        "antsApplyTransforms -d 3 -n MultiLabel -i {input.ribbon} -r {input.ref} -o {output.ribbon} >> {log} 2>&1"

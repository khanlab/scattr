# Directories
freesurfer_dir = str(Path(config["output_dir"]) / "freesurfer")
if config.get("freesurfer_dir"):
    freesurfer_dir = config.get("freesurfer_dir")

log_dir = str(Path(config["output_dir"]) / ".logs" / "freesurfer")

# Licenses
if config.get("fs_license"):
    fs_license = config["fs_license"]
else:
    fs_license = os.getenv("FS_LICENSE")
# Raise error if still no fs_license found
if not fs_license:
    raise ValueError("No Freesurfer license available...exiting")

# BIDS partials
bids_fs_out = partial(
    bids,
    root=freesurfer_dir,
    datatype="anat",
    **inputs.subj_wildcards,
)

bids_log = partial(
    bids,
    root=log_dir,
    **inputs.subj_wildcards,
)

"""Freesurfer references (with additional in rules as necessary)
B. Fischl, A. van der Kouwe, C. Destrieux, E. Halgren, F. Ségonne, D.H. Salat, 
E. Busa, L.J. Seidman, J. Goldstein, D. Kennedy, V. Caviness, N. Makris, 
B. Rosen, A.M. Dale. Automatically parcellating the human cerebral cortex. 
Cereb. Cortex, 14 (2004), pp. 11-22, 10.1093/cercor/bhg087
"""


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
    threads: 1
    resources:
        mem_mb=4000,
        time=10,
    group:
        "dseg_tsv"
    shell:
        "cp -v {input.fs_tsv} {output.fs_tsv}"


rule thalamic_segmentation:
    """Perform thalamus segmentation

    Reference: J.E. Iglesias, R. Insausti, G. Lerma-Usabiaga, M. Bocchetta, 
    K. Van Leemput, D.N. Greve, A. van der Kouwe, B. Fischl, 
    C. Caballero-Gaudes, P.M. Paz-Alonso. A probabilistic atlas of the 
    human thalamic nuclei combining ex vivo MRI and histology. 
    NeuroImage, 183 (2018), pp. 314-326, 10.1016/j.neuroimage.2018.08.012

    NOTE: Output file is defined by Freesurfer script
    """
    input:
        freesurfer_dir=freesurfer_dir,
    params:
        fs_license=fs_license,
    output:
        thal_seg=str(
            Path(bids(root=freesurfer_dir, **inputs.subj_wildcards)).parent
            / "mri"
            / "ThalamicNuclei.v12.T1.mgz"
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        bids_log(suffix="thalamicSegmentation.log"),
    group:
        "freesurfer"
    container:
        config["singularity"]["scattr"]
    shell:
        """
        FS_LICENSE={params.fs_license}
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} 
        export SUBJECTS_DIR={input.freesurfer_dir} 

        mkdir -p {input.freesurfer_dir}/sub-{wildcards.subject}/scripts 

        segmentThalamicNuclei.sh sub-{wildcards.subject} &> {log}
        """


rule mgz2nii:
    """
    Convert from .mgz to .nii.gz

    NOTE: During conversion, files are renamed to BIDS-esque formatting
    """
    input:
        thal=rules.thalamic_segmentation.output.thal_seg
        if not config.get("skip_thal_seg")
        else [],
        aparcaseg=str(
            Path(bids(root=freesurfer_dir, **inputs.subj_wildcards)).parent
            / "mri"
            / "aparc+aseg.mgz"
        ),
    params:
        freesurfer_dir=freesurfer_dir,
        fs_license=fs_license,
    output:
        thal=bids_fs_out(
            space="Freesurfer",
            desc="FreesurferThal",
            suffix="dseg.nii.gz",
        ),
        aparcaseg=bids_fs_out(
            space="Freesurfer",
            desc="aparcaseg",
            suffix="dseg.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
    log:
        bids_log(suffix="mgz2nii.log"),
    group:
        "freesurfer"
    container:
        config["singularity"]["scattr"]
    shell:
        """
        FS_LICENSE={params.fs_license} 
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} 
        export SUBJECTS_DIR={params.freesurfer_dir} 

        mri_convert {input.thal} {output.thal} &> {log} 

        mri_convert {input.aparcaseg} {output.aparcaseg} >> {log} 2>&1 
        """


rule fs_xfm_to_native:
    """Transform from Freesurfer space to T1w space"""
    input:
        thal=rules.mgz2nii.output.thal,
        aparcaseg=rules.mgz2nii.output.aparcaseg,
        ref=lambda wildcards: expand(
            inputs["T1w"].path,
            zip,
            **filter_list(inputs["T1w"].input_zip_lists, wildcards)
        )[0],
    output:
        thal=bids_fs_out(
            space="T1w",
            desc="FreesurferThal",
            suffix="dseg.nii.gz",
        ),
        aparcaseg=bids_fs_out(
            space="T1w",
            desc="aparcaseg",
            suffix="dseg.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        bids_log(suffix="fsXfmToNative.log"),
    group:
        "freesurfer"
    container:
        config["singularity"]["scattr"]
    shell:
        """
        antsApplyTransforms -d 3 -n MultiLabel \\
            -i {input.thal} -r {input.ref} -o {output.thal} &> {log} 

        antsApplyTransforms -d 3 -n MultiLabel \\
            -i {input.aparcaseg} -r {input.ref} \\
            -o {output.aparcaseg} >> {log} 2>&1"
        """

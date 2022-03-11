from os.path import join 

rule thalamic_segmentation:
    input:
        freesurfer_dir=join(config['bids_dir'], 'derivatives/freesurfer'),
    params:
        fs_license=config["fs_license"],
        threads=workflow.cores,
    output:
        thal_seg=bids(
            root=join(config["output_dir"], "freesurfer"),
            datatype="anat",
            suffix="T1.mgz"
            space="Freesurfer",
            **config["subj_wildcards"],
        ),
    container: 
        config['singularity']['freesurfer']
    shell:
        "FS_LICENSE={params.fs_license} "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={params.threads} "
        "SUBJECTS_DIR={input.freesurfer_dir} "
        "segmentThalamicNuclei.sh {{subject}}"

rule mgz_to_nii:
    input: 
        thal=rules.thalamic_segmentation.output.thal_seg,
        aparcaseg=join(config["output_dir"], "freesurfer/{subject}/mri/aparc+aseg.mgz")
    params:
        fs_license=config["fs_license"],
        threads=workflow.cores,
    output:
        thal=bids(
            root=join(config["output_dir"], "freesurfer"),
            datatype="anat",
            suffix="T1.nii.gz"
            space="Freesurfer",
            **config["subj_wildcards"],
        ),
        aparcaseg=bids(
            root=join(config["output_dir"], "freesurfer"),
            datatype="anat",
            suffix="aparcaseg.nii.gz"
            space="Freesurfer",
            **config["subj_wildcards"],
        )
    container:
        config['singularity']['freesurfer']
    shell:
        "FS_LICENSE={params.fs_license} "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={params.threads} "
        "mri_convert {input.thal} {output.thal} "
        "mri_convert {input.aparcaseg} {output.aparcaseg} "

rule fs_xfm_to_native:
    input:
        thal=rules.mgz_to_nii.output.thal,
        aparcaseg=rules.mgz_to_nii.output.aparcaseg,
        ref=bids(
            root=config["bids_dir"],
            datatype="anat",
            suffix="T1w.nii.gz",
            **config["subj_wildcards"],
        ),
    output:
        thal=bids(
            root=join(config["output_dir"], "freesurfer"),
            datatype="anat",
            suffix="T1.nii.gz"
            space="T1w",
            **config["subj_wildcards"],
        ),
        aparcaseg=bids(
            root=join(config["output_dir"], 'freesurfer'),
            datatype='anat',
            suffix="aparcaseg.nii.gz",
            space='T1w',
            **config["subj_wildcards"],
        )
    container:
        config['singularity']['neuroglia-core'],
    shell:
        "antsApplyTransforms -d 3 -n MultiLabel -i {input.thal} -r {input.ref} -o {output.thal} "
        "antsApplyTransforms -d 3 -n MultiLabel -i {input.aparcaseg} -r {input.ref} -o {output.aparcaseg} "
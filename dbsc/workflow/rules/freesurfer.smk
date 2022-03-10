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
        "segmentThalamicNuclei.sh {{subject:4}}"

rule mgz_to_nii:
    input: 
        mgz=rule.thalamic_segmentation.output.thal_seg,
    params:
        fs_license=config["fs_license"],
        threads=workflow.cores,
    output:
        nii=bids(
            root=join(config["output_dir"], "freesurfer"),
            datatype="anat",
            suffix="T1.nii.gz"
            space="Freesurfer",
            **config["subj_wildcards"],
        ),
    container:
        config['singularity']['freesurfer']
    shell:
        "FS_LICENSE={params.fs_license} "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={params.threads} "
        "mri_convert {input.mgz} {output.nii}"

rule xfm_to_native:
    input:
        nii=rule.mgz_to_nii.output.nii,
        ref=bids(
            root=config["bids_dir"],
            datatype="anat",
            suffix="T1w.nii.gz",
            **config["subj_wildcards"],
        ),
    output:
        nii=bids(
            root=join(config["output_dir"], "freesurfer"),
            datatype="anat",
            suffix="T1.nii.gz"
            space="Native",
            **config["subj_wildcards"],
        ),
    container:
        config['singularity']['neuroglia-core'],
    shell:
        "antsApplyTransforms -d 3 -n MultiLabel -i {input.nii} -r {input.ref} -o {output.nii}"
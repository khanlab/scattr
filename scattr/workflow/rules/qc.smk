bids_labelmerge = partial(
    bids,
    root=str(Path(labelmerge_dir) / "combined"),
    **config["subj_wildcards"],
)

rule overlay_image:
    input:
        qc_labels = bids_labelmerge(
                space="T1w",
                desc="combined",
                suffix="dseg.nii.gz",
            ), 
        t1w_image = config["input_path"]["T1w"]

    output:
        segment_overlay = bids(root =  str(Path(config["output_dir"]) / "qc_dir"), 
        datatype = 'anat',
        space = 'T1w',
        suffix='qc_img.png',
        **config["subj_wildcards"]
        )

    script:
        "scripts/qc/segmentation_image.py" 
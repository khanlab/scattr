# BIDS partials
bids_labelmerge = partial(
    bids,
    root=str(Path(labelmerge_dir) / "combined"),
    **config["subj_wildcards"],
)

bids_qc = partial(
    bids,
    root=str(Path(config["output_dir"] / "qc")),
    datatype="anat",
    space="T1w",
    **config["subj_wildcards"],
)


rule overlay_segment:
    input:
        qc_labels=bids_labelmerge(
            space="T1w",
            desc="combined",
            suffix="dseg.nii.gz",
        ),
        t1w_image=config["input_path"]["T1w"],
    output:
        qc_png=bids_qc(
            desc="labelmerge",
            suffix="dseg.png",
        ),
        qc_html=bids_qc(
            desc="labelmerge",
            suffix="dseg.html",
        ),
    script:
        "../scripts/qc/segmentation_image.py"

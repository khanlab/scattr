# BIDS partials
bids_labelmerge = partial(
    bids,
    root=str(Path(labelmerge_dir) / "combined"),
    **config["subj_wildcards"],
)

bids_qc = partial(
    bids,
    root=str(Path(config["output_dir"]) / "qc"),
    datatype="anat",
    space="T1w",
    **config["subj_wildcards"],
)


rule segment_qc:
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
            suffix="dsegQC.png",
        ),
        qc_html=bids_qc(
            desc="labelmerge",
            suffix="dsegQC.html",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=15,
    group:
        "qc"
    container:
        config["singularity"]["scattr"]
    script:
        "../scripts/qc/segmentation_qc.py"


rule registration_qc:
    input:
        template_t1w=rules.reg2native.output.t1w_nativespace,
        native_t1w=config["input_path"]["T1w"],
    output:
        qc_png=bids_qc(
            desc=f"{config['Space']}toNative",
            suffix="regQC.png",
        ),
        qc_html=bids_qc(
            desc=f"{config['Space']}toNative",
            suffix="regQC.html",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=15,
    group:
        "qc"
    container:
        config["singularity"]["scattr"]
    script:
        "../scripts/qc/registration_qc.py"


rule gather_qc:
    input:
        dseg_png=expand(
            rules.segment_qc.output.qc_png,
            subject=config["input_lists"]["T1w"]["subject"],
        ),
        dseg_html=expand(
            rules.segment_qc.output.qc_html,
            subject=config["input_lists"]["T1w"]["subject"],
        ),
        reg_png=expand(
            rules.registration_qc.output.qc_png,
            subject=config["input_lists"]["T1w"]["subject"],
        ),
        reg_html=expand(
            rules.registration_qc.output.qc_html,
            subject=config["input_lists"]["T1w"]["subject"],
        ),

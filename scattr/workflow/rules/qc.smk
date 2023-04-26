# BIDS partials
bids_labelmerge = partial(
    bids,
    root=str(Path(labelmerge_dir) / "combined")
    if not config.get("skip_labelmerge")
    else config.get("labelmerge_base_dir") or zona_dir,
    **inputs.subj_wildcards,
)

bids_qc = partial(
    bids,
    root=str(Path(config["output_dir"]) / "qc"),
    datatype="anat",
    space="T1w",
    **inputs.subj_wildcards,
)


rule segment_qc:
    input:
        qc_labels=bids_labelmerge(
            space="T1w",
            datatype=None if not config.get("skip_labelmerge") else "anat",
            desc="combined"
            if not config.get("skip_labelmerge")
            else config.get("labelmerge_base_desc"),
            suffix="dseg.nii.gz",
        ),
        t1w_image=inputs["T1w"].path,
    output:
        qc_png=bids_qc(
            desc="labelmerge"
            if not config.get("skip_labelmerge")
            else config.get("labelmerge_base_desc"),
            suffix="dsegQC.png",
        ),
        qc_html=bids_qc(
            desc="labelmerge"
            if not config.get("skip_labelmerge")
            else config.get("labelmerge_base_desc"),
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
        moving_nii=rules.reg2native.output.t1w_nativespace,
        fixed_nii=inputs["T1w"].path,
    params:
        cuts=7,
    output:
        qc_svg=bids_qc(
            desc=f"{config['Space']}toNative",
            suffix="regQC.svg",
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
            zip,
            **inputs["T1w"].input_zip_lists,
        ),
        dseg_html=expand(
            rules.segment_qc.output.qc_html,
            zip,
            **inputs["T1w"].input_zip_lists,
        ),
        reg_png=expand(
            rules.registration_qc.output.qc_svg,
            zip,
            **inputs["T1w"].input_zip_lists,
        ),
        reg_html=expand(
            rules.registration_qc.output.qc_html,
            zip,
            **inputs["T1w"].input_zip_lists,
        ),

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
        # HOW TO PASS APPROPRIATE T1W?
        t1w_image=lambda wildcards: expand(
            inputs["T1w"].path,
            zip,
            **filter_list(inputs["T1w"].input_zip_lists, wildcards)
        )[0],
    output:
        qc_png=report(
            bids_qc(
                desc="labelmerge"
                if not config.get("skip_labelmerge")
                else config.get("labelmerge_base_desc"),
                suffix="dsegQC.png",
            ),
            caption="../report/dseg_qc.rst",
            category="Label Segmentation QC",
            subcategory="Static",
        ),
        qc_html=report(
            bids_qc(
                desc="labelmerge"
                if not config.get("skip_labelmerge")
                else config.get("labelmerge_base_desc"),
                suffix="dsegQC.html",
            ),
            category="Label Segmentation QC",
            subcategory="Interactive",
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
        fixed_nii=lambda wildcards: expand(
            inputs["T1w"].path,
            zip,
            **filter_list(inputs["T1w"].input_zip_lists, wildcards)
        )[0],
    params:
        cuts=7,
    output:
        qc_svg=report(
            bids_qc(
                desc=f"{config['Space']}toNative",
                suffix="regQC.svg",
            ),
            caption="../report/reg_qc.rst",
            category="Template to Subject Registration",
            subcategory="Static GIF",
        ),
        qc_html=report(
            bids_qc(
                desc=f"{config['Space']}toNative",
                suffix="regQC.html",
            ),
            category="Template to Subject Registration",
            subcategory="Interactive",
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
            **inputs["T1w"].input_zip_lists
        ),
        dseg_html=expand(
            rules.segment_qc.output.qc_html,
            zip,
            **inputs["T1w"].input_zip_lists
        ),
        reg_svg=expand(
            rules.registration_qc.output.qc_svg,
            zip,
            **inputs["T1w"].input_zip_lists
        ),
        reg_html=expand(
            rules.registration_qc.output.qc_html,
            zip,
            **inputs["T1w"].input_zip_lists
        ),

# BIDS partials
bids_labelmerge = partial(
    bids,
    root=str(Path(labelmerge_dir) / "combined")
    if not config.get("skip_labelmerge")
    else config.get("labelmerge_base_dir") or zona_dir,
    **inputs_t1w.subj_wildcards,
)

bids_qc = partial(
    bids,
    root=str(Path(config["output_dir"]) / "qc"),
    datatype="anat",
    space="T1w",
    **inputs_t1w.subj_wildcards,
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
        t1w_image=lambda wildcards: inputs_t1w["T1w"]
        .filter(**wildcards)
        .expand()[0],
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
    conda:
        "../envs/neurovis.yaml"
    script:
        "../scripts/qc/segmentation_qc.py"


rule registration_qc:
    input:
        moving_nii=rules.reg2native.output.t1w_nativespace,
        fixed_nii=lambda wildcards: inputs_t1w["T1w"]
        .filter(**wildcards)
        .expand()[0],
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
    conda:
        "../envs/neurovis.yaml"
    script:
        "../scripts/qc/registration_qc.py"


rule gather_qc:
    input:
        dseg_png=inputs_t1w["T1w"].expand(rules.segment_qc.output.qc_png),
        dseg_html=inputs_t1w["T1w"].expand(rules.segment_qc.output.qc_html),
        reg_svg=inputs_t1w["T1w"].expand(rules.registration_qc.output.qc_svg),
        reg_html=inputs_t1w["T1w"].expand(rules.registration_qc.output.qc_html),

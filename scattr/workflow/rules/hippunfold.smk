from pathlib import Path
import glob
import nibabel and nib

hippunfold_root = Path(config["bids_dir"]) / "derivatives" / "hippunfold"

hippunfold_overlay_dir = (
    config.get("hippunfold_overlay_dir")
    or str(Path(config["output_dir"]) / "hippunfold_bids")
)

log_dir = Path(config["output_dir"]) / ".logs" / "hippunfold"
Path(hippunfold_overlay_dir).mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

bids_hippu_out = partial(
	bids,
	root=hippunfold_overlay_dir,
	datatype='anat',
	**inputs_t1w.subj_wildcards,
)

bids_hippu_log = partial(
    bids,
    root=str(log_dir),
    **inputs_t1w.subj_wildcards,
)


rule cp_hippu_tsv:
    """
    Copy HippUnfold TSV (provided by user in config)
    into hippunfold overlay directory.
    """
    input:
        hip_tsv=str(
            Path(workflow.basedir).parent
            / Path(config["hippunfold"]["tsv"])
        )
    output:
        hip_tsv=str(
            Path(hippunfold_overlay_dir)
            / "desc-HippUnfoldSubfields_dseg.tsv"
        )
    threads: 1
    resources:
        mem_mb=1000
    shell:
        "cp -v {input.hip_tsv} {output.hip_tsv}"


rule hippunfold_merge_subfields:
    """
    Locate the HippUnfold left/right dseg outputs from:
    <bids_dir>/derivatives/hippunfold/sub-XXX/anat

    Merge into a single:
    <output_dir>/hippunfold_bids/sub-XXX/anat/...desc-HippUnfoldSubfields_dseg.nii.gz
    """
    input:
        t1=inputs_t1w["T1w"]
    output:
        dseg=bids_hippu_out(
            space="T1w",
            desc="HippUnfoldSubfields",
            suffix="dseg.nii.gz",
        )
    threads: 2
    resources:
        mem_mb=8000
    log:
        bids_hippu_log(suffix="merge.log")
    run:
        import nibabel as nib
        import glob
        from pathlib import Path

        sub = wildcards.subject

        base = hippunfold_root / f"sub-{sub}" / "anat"

        patt_L = base / f"sub-{sub}_hemi-L_space-T1w_label-hipp_desc-subfields_atlas-multihist7_dseg.nii.gz"
        patt_R = base / f"sub-{sub}_hemi-R_space-T1w_label-hipp_desc-subfields_atlas-multihist7_dseg.nii.gz"

        hits_L = glob.glob(str(patt_L))
        hits_R = glob.glob(str(patt_R))

        if not hits_L or not hits_R:
            raise FileNotFoundError(
                f"Missing HippUnfold outputs:\n  {patt_L}\n  {patt_R}"
            )

        imgL = nib.load(hits_L[0])
        imgR = nib.load(hits_R[0])
        L = imgL.get_fdata().astype(int)
        R = imgR.get_fdata().astype(int)

        R[R > 0] += 100

        merged = L.copy()
        nz = R > 0
        merged[nz] = R[nz]

        Path(str(output.dseg)).parent.mkdir(parents=True, exist_ok=True)
        nib.Nifti1Image(merged, imgL.affine, imgL.header).to_filename(str(output.dseg))

        with open(str(log), "w") as LOG:
            LOG.write(
                f"Merged:\nL={hits_L[0]}\nR={hits_R[0]}\n→ {output.dseg}\n"
            )
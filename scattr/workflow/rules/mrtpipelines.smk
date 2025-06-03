import nibabel as nib
import numpy as np


# Directories
responsemean_dir = config.get("responsemean_dir")
mrtrix_dir = str(Path(config["output_dir"]) / "mrtrix")
labelmerge_dir = str(Path(config["output_dir"]) / "labelmerge")
zona_dir = str(Path(config["output_dir"]) / "zona_bb_subcortex")
log_dir = str(Path(config["output_dir"]) / ".logs" / "mrtrix")

# Make directory if it doesn't exist
Path(mrtrix_dir).mkdir(parents=True, exist_ok=True)

# Parameters
shells = config.get("shells")
lmax = config.get("lmax")


# BIDS partials
bids_response_out = partial(
    bids,
    root=mrtrix_dir,
    datatype="response",
    suffix="response.txt",
)

bids_dti_out = partial(
    bids,
    root=mrtrix_dir,
    datatype="dti",
    model="dti",
    **inputs_dwi.subj_wildcards,
)

bids_tractography_out = partial(
    bids,
    root=mrtrix_dir,
    datatype="tractography",
    **inputs_dwi.subj_wildcards,
)

bids_anat_out = partial(
    bids,
    root=mrtrix_dir,
    datatype="anat",
    **inputs_t1w.subj_wildcards,
)

bids_labelmerge = partial(
    bids,
    root=str(Path(labelmerge_dir) / "combined")
    if not config.get("skip_labelmerge")
    else config.get("labelmerge_base_dir") or zona_dir,
    **inputs_t1w.subj_wildcards,
)

bids_log = partial(
    bids,
    root=log_dir,
    **inputs_dwi.subj_wildcards,
)

"""Mrtrix3 reference (additional citations are included per rule as necessary):
Tournier, J.-D.; Smith, R. E.; Raffelt, D.; Tabbara, R.; Dhollander, T.; 
Pietsch, M.; Christiaens, D.; Jeurissen, B.; Yeh, C.-H. & Connelly, A. 
MRtrix3: A fast, flexible and open software framework for medical image 
processing and visualisation. NeuroImage, 2019, 202, 116137
"""


include: "mrtpipelines/preproc.smk"
include: "mrtpipelines/tractography.smk"
include: "mrtpipelines/filter_tck.smk"


rule dwi2tensor:
    input:
        dwi=rules.dwinormalise.output.dwi,
        mask=rules.nii2mif.output.mask,
    output:
        dti=bids_dti_out(suffix="tensor.mif"),
        fa=bids_dti_out(suffix="fa.mif"),
        ad=bids_dti_out(suffix="ad.mif"),
        rd=bids_dti_out(suffix="rd.mif"),
        md=bids_dti_out(suffix="md.mif"),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        bids_log(suffix="dwi2tensor.log"),
    group:
        "dwiproc"
    container:
        config["singularity"]["scattr"]
    conda:
        "../envs/mrtrix3.yaml"
    shell:
        """
        dwi2tensor -nthreads {threads} -mask {input.mask} \\
        {input.dwi} {output.dti} &> {log}

        tensor2metric -nthreads {threads} -mask {input.mask} {output.dti} \\
        -fa {output.fa} -ad {output.ad} -rd {output.rd} \\
        -adc {output.md} >> {log} 2>&1
        """

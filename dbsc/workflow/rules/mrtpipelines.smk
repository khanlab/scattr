from pathlib import Path
from functools import partial

import numpy as np


# Directories
mrtrix_dir = str(Path(config["output_dir"]) / "mrtrix")

# Make directory if it doesn't exist
Path(mrtrix_dir).mkdir(parents=True, exist_ok=True)

# Parameters
responsemean_flag = config.get("responsemean_dir")
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
    **config["subj_wildcards"],
)

bids_tractography_out = partial(
    bids,
    root=mrtrix_dir,
    datatype="tractography",
    **config["subj_wildcards"],
)

bids_anat_out = partial(
    bids,
    root=mrtrix_dir,
    datatype="anat",
    **config["subj_wildcards"],
)

# Mrtrix3 citation (additional citations are included per rule as necessary):
# Tournier, J.-D.; Smith, R. E.; Raffelt, D.; Tabbara, R.; Dhollander, T.; Pietsch, M.; Christiaens, D.; Jeurissen, B.; Yeh, C.-H. & Connelly, A. MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. NeuroImage, 2019, 202, 116137


# ------------ MRTRIX PREPROC BEGIN ----------#
rule nii2mif:
    input:
        dwi=config["input_path"]["dwi"],
        bval=re.sub(".nii.gz", ".bval", config["input_path"]["dwi"]),
        bvec=re.sub(".nii.gz", ".bvec", config["input_path"]["dwi"]),
        mask=config["input_path"]["mask"],
    output:
        dwi=bids(
            root=mrtrix_dir,
            datatype="dwi",
            suffix="dwi.mif",
            **config["subj_wildcards"]
        ),
        mask=bids(
            root=mrtrix_dir,
            datatype="dwi",
            suffix="brainmask.mif",
            **config["subj_wildcards"]
        ),
    threads: workflow.cores
    group:
        "subject_1"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrconvert -nthreads {threads} -fslgrad {input.bvec} {input.bval} {input.dwi} {output.dwi} && "
        "mrconvert -nthreads {threads} {input.mask} {output.mask}"


rule dwi2response:
    """
    Estimate response functions using Dhollander algorithm

    Dhollander, T.; Mito, R.; Raffelt, D. & Connelly, A. Improved white matter response function estimation for 3-tissue constrained spherical deconvolution. Proc Intl Soc Mag Reson Med, 2019, 555
    """
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
    params:
        shells=f"-shells {shells}" if shells else "",
        lmax=f"-lmax {lmax}" if lmax else "",
    output:
        wm_rf=bids_response_out(
            desc="wm",
            **config["subj_wildcards"],
        ),
        gm_rf=bids_response_out(
            desc="gm",
            **config["subj_wildcards"],
        ),
        csf_rf=bids_response_out(
            desc="csf",
            **config["subj_wildcards"],
        ),
    threads: workflow.cores
    group:
        "subject_1"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2response dhollander {input.dwi} {output.wm_rf} {output.gm_rf} {output.csf_rf} -nthreads {threads} -mask {input.mask} {params.shells} {params.lmax} "


rule responsemean:
    """Compute average response function"""
    input:
        subject_rf=bids_response_out(
            desc="{tissue}",
            **config["subj_wildcards"],
        ),
    output:
        avg_rf=bids_response_out(
            root=str(Path(mrtrix_dir) / "avg"),
            desc="{tissue}",
            **config["subj_wildcards"],
        ),
    threads: workflow.cores
    group:
        "group"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "responsemean {input.subject_rf} {output.avg_rf} -nthreads {threads}"


rule dwi2fod:
    """Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data. NeuroImage, 2014, 103, 411-426"""
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
        wm_rf=(
            str(Path(config["responsemean_dir"]) / "desc-wm_response.txt")
            if responsemean_flag
            else expand(
                rules.responsemean.output.avg_rf, tissue="wm", allow_missing=True
            )
        ),
        gm_rf=(
            str(Path(config["responsemean_dir"]) / "desc-gm_response.txt")
            if responsemean_flag
            else expand(
                rules.responsemean.output.avg_rf, tissue="gm", allow_missing=True
            )
        ),
        csf_rf=(
            str(Path(config["responsemean_dir"]) / "desc-csf_response.txt")
            if responsemean_flag
            else expand(
                rules.responsemean.output.avg_rf, tissue="csf", allow_missing=True
            )
        ),
    params:
        shells=f"-shells {shells}" if shells else "",
    output:
        wm_fod=bids_response_out(
            model="csd",
            desc="wm",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
        gm_fod=bids_response_out(
            model="csd",
            desc="gm",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
        csf_fod=bids_response_out(
            model="csd",
            desc="csf",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
    threads: workflow.cores
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2fod -nthreads {threads} {params.shells} -mask {input.mask} msmt_csd {input.dwi} {input.wm_rf} {output.wm_fod} {input.gm_rf} {output.gm_fod} {input.csf_rf} {output.csf_fod}"


rule mtnormalise:
    """
    Raffelt, D.; Dhollander, T.; Tournier, J.-D.; Tabbara, R.; Smith, R. E.; Pierre, E. & Connelly, A. Bias Field Correction and Intensity Normalisation for Quantitative Analysis of Apparent Fibre Density. In Proc. ISMRM, 2017, 26, 3541

    Dhollander, T.; Tabbara, R.; Rosnarho-Tornstrand, J.; Tournier, J.-D.; Raffelt, D. & Connelly, A. Multi-tissue log-domain intensity and inhomogeneity normalisation for quantitative apparent fibre density. In Proc. ISMRM, 2021, 29, 2472
    """
    input:
        wm_fod=rules.dwi2fod.output.wm_fod,
        gm_fod=rules.dwi2fod.output.gm_fod,
        csf_fod=rules.dwi2fod.output.csf_fod,
        mask=rules.nii2mif.output.mask,
    output:
        wm_fod=bids_response_out(
            model="csd",
            desc="wm",
            suffix="fodNormalized.mif",
            **config["subj_wildcards"],
        ),
        gm_fod=bids_response_out(
            model="csd",
            desc="gm",
            suffix="fodNormalized.mif",
            **config["subj_wildcards"],
        ),
        csf_fod=bids_response_out(
            model="csd",
            desc="csf",
            suffix="fodNormalized.mif",
            **config["subj_wildcards"],
        ),
    threads: workflow.cores
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mtnormalise -nthreads {threads} -mask {input.mask} {input.wm_fod} {output.wm_fod} {input.gm_fod} {output.gm_fod} {input.csf_fod} {output.csf_fod}"


# DTI (Tensor) Processing
rule dwinormalise:
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
    output:
        dwi=bids(
            root=mrtrix_dir,
            datatype="dwi",
            desc="normalized",
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),
    threads: workflow.cores
    container:
        config["singularity"]["mrtrix"]
    group:
        "subject_1"
    shell:
        "dwinormalise individual -nthreads {threads} {input.dwi} {input.mask} {output.dwi}"


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
    threads: workflow.cores
    group:
        "subject_1"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2tensor -nthreads {threads} -mask {input.mask} {input.dwi} {output.dti} && "
        "tensor2metric -nthreads {threads} -mask {input.mask} {output.dti} -fa {output.fa} -ad {output.ad} -rd {output.rd} -adc {output.md}"


# --------------- MRTRIX PREPROC END --------------#


# ------------ MRTRIX TRACTOGRAPHY BEGIN ----------#
rule tckgen:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        fod=rules.mtnormalise.output.wm_fod,
        mask=rules.nii2mif.output.mask,
        cortical_ribbon=rules.fs_xfm_to_native.output.ribbon,
        convex_hull=rules.create_convex_hull.output.convex_hull,
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
    params:
        step=config["step"],
        sl=config["sl_count"],
    output:
        tck=bids_tractography_out(
            desc="iFOD2",
            suffix="tractography.tck",
        ),
    threads: 32
    resources:
        mem_mb=128000,
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tckgen -nthreads {threads} -algorithm iFOD2 -step {params.step} -select {params.sl} -exclude {input.cortical_ribbon} -exclude {input.convex_hull} -include {input.subcortical_seg} -mask {input.mask} -seed_image {input.mask} {input.fod} {output.tck}"


rule tcksift2:
    # Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The effects of SIFT on the reproducibility and biological accuracy of the structural connectome. NeuroImage, 2015, 104, 253-265
    input:
        tck=rules.tckgen.output.tck,
        fod=rules.mtnormalise.output.wm_fod,
    output:
        weights=bids_tractography_out(
            desc="iFOD2",
            suffix="tckWeights.txt",
        ),
        mu=bids_tractography_out(desc="iFOD2", suffix="muCoefficient.txt"),
    threads: 32
    resources:
        mem_mb=128000,
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tcksift2 -nthreads {threads} -out_mu {output.mu} {input.tck} {input.fod} {output.weights}"


# TODO (v0.2): ADD OPTION TO OUTPUT TDI MAP


rule tck2connectome:
    """
    Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The effects of SIFT on the reproducibility and biological accuracy of the structural connectome. NeuroImage, 2015, 104, 253-265"

    TODO (v0.2): ADD OPTION TO MULTIPLY BY MU COEFFICIENT
    """
    input:
        weights=rules.tcksift2.output.weights,
        tck=rules.tckgen.output.tck,
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
    params:
        radius=config["radial_search"],
    output:
        sl_assignment=bids_tractography_out(
            desc="subcortical",
            suffix="nodeAssignment.txt",
        ),
        node_weights=bids_tractography_out(
            desc="subcortical",
            suffix="nodeWeights.csv",
        ),
    threads: 32
    resources:
        mem_mb=128000,
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tck2connectome -nthreads {threads} -zero_diagonal -stat_edge sum -assignment_radial_search {params.radius} -tck_weights_in {input.weights} -out_assignments {output.sl_assignment} -symmetric {input.tck} {input.subcortical_seg} {output.node_weights} "


rule connectome2tck:
    input:
        node_weights=rules.tck2connectome.output.node_weights,
        sl_assignment=rules.tck2connectome.output.sl_assignment,
        tck=rules.tckgen.output.tck,
    params:
        nodes=",".join(str(num) for num in range(1, 73)),
    output:
        edge_weight=temp(
            bids_tractography_out(
                desc="subcortical",
                suffix="tckWeights",
            )
        ),
        edge_tck=temp(
            bids_tractography_out(
                desc="subcortical",
                suffix="from",
            )
        ),
    threads: 32
    resources:
        mem_mb=128000,
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "connectome2tck -nthreads {threads} -nodes {params.nodes} -exclusive -filters_per_edge -tck_weights_in {input.node_weights} -prefix_tck_weights_out {output.edge_weight} {input.tck} {input.sl_assignment} {output.edge_tck} "


# NOTE: Use labelmerge split segs here?
rule create_roi_mask:
    # EXPAND OVER NODE1 IN RULE ALL
    input:
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
    output:
        roi_mask=temp(
            bids_anat_out(
                desc="{node1}",
                suffix="mask.mif",
            )
        ),
    threads: workflow.cores
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrcalc -nthreads {threads} {input.subcortical_seg} {wildcards.node1} -eq {output.roi_mask}"


rule create_exclude_mask:
    input:
        roi1=bids_anat_out(
            desc="{node1}",
            suffix="mask.mif",
        ),
        roi2=bids_anat_out(
            desc="{node2}",
            suffix="mask.mif",
        ),
        lZI=bids_anat_out(desc="21", suffix="mask.mif"),
        rZI=bids_anat_out(desc="22", suffix="mask.mif"),
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
    output:
        filter_mask=temp(
            bids_anat_out(
                desc="exclude{node1}AND{node2}",
                suffix="mask.mif",
            )
        ),
    threads: workflow.cores
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrcalc -nthreads {threads} {input.subcortical_seg} 0 -neq {input.roi1} -sub {input.roi2} -sub {input.lZI} -sub {input.rZI} -sub {output.filter_mask}"


rule filter_tck:
    input:
        filter_mask=rules.create_exclude_mask.output.filter_mask,
        tck=rules.connectome2tck.output.edge_tck,
        weights=rules.tck2connectome.output.node_weights,
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
    output:
        filtered_tck=temp(
            bids_tractography_out(
                desc="from{node1}to{node2}",
                suffix="tractography.tck",
            ),
        ),
        filtered_weights=temp(
            bids_tractography_out(
                desc="from{node1}to{node2}",
                suffix="weights.csv",
            ),
        ),
    threads: 32
    resources:
        mem_mb=128000,
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tckedit -nthreads {threads} -exclude {input.filter_mask} -tck_weights_in {input.weights} -tck_weights_out {output.filtered_weights} {input.tck} {output.filtered_tck}"


idxes = np.triu_indices(72, k=1)


rule combine_filtered:
    input:
        tck=expand(
            rules.filter_tck.output.filtered_tck,
            zip,
            node1=list(idxes[0] + 1),
            node2=list(idxes[1] + 1),
            allow_missing=True,
        ),
        weights=expand(
            rules.filter_tck.output.filtered_weights,
            zip,
            node1=list(idxes[0] + 1),
            node2=list(idxes[1] + 1),
            allow_missing=True,
        ),
    output:
        combined_tck=bids_tractography_out(
            desc="filteredsubcortical",
            suffix="tractography.tck",
        ),
        combined_weights=bids_tractography_out(
            desc="filteredsubcortical",
            suffix="tckWeights.txt",
        ),
    threads: workflow.cores
    resources:
        mem_mb=128000,
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tckedit {input.tck} {output.combined_tck} && "
        "cat {input.weights} >> {output.combined_weights} "


rule filtered_tck2connectome:
    input:
        weights=rules.combine_filtered.output.combined_weights,
        tck=rules.combine_filtered.output.combined_tck,
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
    params:
        radius=config["radial_search"],
    output:
        sl_assignment=bids_tractography_out(
            desc="filteredsubcortical",
            suffix="nodeAssignment.txt",
        ),
        node_weights=bids_tractography_out(
            desc="filteredsubcortical",
            suffix="nodeWeights.csv",
        ),
    threads: 32
    resources:
        mem_mb=128000,
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tck2connectome -nthreads {threads} -zero_diagonal -stat_edge sum -assignment_radial_search {params.radius} -tck_weights_in {input.weights} -out_assignments {output.sl_assignment} -symmetric {input.tck} {input.subcortical_seg} {output.node_weights} -force"


# ------------ MRTRIX TRACTOGRAPHY END ----------#

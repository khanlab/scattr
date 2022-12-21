import nibabel as nib 
import numpy as np


# Directories
responsemean_dir = config.get("responsemean_dir")
dwi_dir = config.get("dwi_dir")
mrtrix_dir = str(Path(config["output_dir"]) / "mrtrix")
labelmerge_dir = str(Path(config["output_dir"]) / "labelmerge")

# Make directory if it doesn't exist
Path(mrtrix_dir).mkdir(parents=True, exist_ok=True)

# Parameters
shells = config.get("shells")
lmax = config.get("lmax")

# BIDS partials
bids_dwi = partial(
    bids,
    root=dwi_dir,
    datatype="dwi",
    space="T1w",
    desc="preproc",
    **config["subj_wildcards"],
)

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

bids_labelmerge = partial(
    bids,
    root=str(Path(labelmerge_dir) / "combined"),
    **config["subj_wildcards"],
)

# Mrtrix3 citation (additional citations are included per rule as necessary):
# Tournier, J.-D.; Smith, R. E.; Raffelt, D.; Tabbara, R.; Dhollander, T.; Pietsch, M.; Christiaens, D.; Jeurissen, B.; Yeh, C.-H. & Connelly, A. MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. NeuroImage, 2019, 202, 116137


# ------------ MRTRIX PREPROC BEGIN ----------#
if dwi_dir:
    print(f"Searching {dwi_dir} for dwi and mask images...")


rule nii2mif:
    input:
        dwi=(
            bids_dwi(suffix="dwi.nii.gz")
            if dwi_dir
            else config["input_path"]["dwi"]
        ),
        bval=(
            bids_dwi(suffix="dwi.bval")
            if dwi_dir
            else re.sub(".nii.gz", ".bval", config["input_path"]["dwi"])
        ),
        bvec=(
            bids_dwi(suffix="dwi.bvec")
            if dwi_dir
            else re.sub(".nii.gz", ".bvec", config["input_path"]["dwi"])
        ),
        mask=(
            bids_dwi(suffix="mask.nii.gz")
            if dwi_dir
            else config["input_path"]["mask"]
        ),
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
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
    log:
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/nii2mif.log",
    group: "dwiproc"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrconvert -nthreads {threads} -fslgrad {input.bvec} {input.bval} {input.dwi} {output.dwi} &> {log} && "
        "mrconvert -nthreads {threads} {input.mask} {output.mask} >> {log} 2>&1"


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
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/dwi2response.log",
    group: "dwiproc"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2response dhollander {input.dwi} {output.wm_rf} {output.gm_rf} {output.csf_rf} -nthreads {threads} -mask {input.mask} {params.shells} {params.lmax} &> {log}"


rule responsemean:
    """Compute average response function"""
    input:
        subject_rf=expand(
            bids_response_out(
                subject="{subject}",
                desc="{tissue}",
            ),
            allow_missing=True,
            subject=config["input_lists"]["T1w"]["subject"],
        ),
    output:
        avg_rf=bids_response_out(
            root=str(Path(mrtrix_dir) / "avg"),
            desc="{tissue}",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
    log:
        f"{config['output_dir']}/logs/mrtrix/{{tissue}}_responsemean.log",
    group: "dwiproc_group"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "responsemean {input.subject_rf} {output.avg_rf} -nthreads {threads} &> {log}"


rule dwi2fod:
    """Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data. NeuroImage, 2014, 103, 411-426"""
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
        wm_rf=(
            str(Path(responsemean_dir) / "desc-wm_response.txt")
            if responsemean_dir
            else expand(
                rules.responsemean.output.avg_rf,
                tissue="wm",
                allow_missing=True,
            )
        ),
        gm_rf=(
            str(Path(responsemean_dir) / "desc-gm_response.txt")
            if responsemean_dir
            else expand(
                rules.responsemean.output.avg_rf,
                tissue="gm",
                allow_missing=True,
            )
        ),
        csf_rf=(
            str(Path(responsemean_dir) / "desc-csf_response.txt")
            if responsemean_dir
            else expand(
                rules.responsemean.output.avg_rf,
                tissue="csf",
                allow_missing=True,
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
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/dwi2fod.log",
    group: "diffmodel"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2fod -nthreads {threads} {params.shells} -mask {input.mask} msmt_csd {input.dwi} {input.wm_rf} {output.wm_fod} {input.gm_rf} {output.gm_fod} {input.csf_rf} {output.csf_fod} &> {log}"


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
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/mtnormalise.log",
    group: "diffmodel"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mtnormalise -nthreads {threads} -mask {input.mask} {input.wm_fod} {output.wm_fod} {input.gm_fod} {output.gm_fod} {input.csf_fod} {output.csf_fod} &> {log}"


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
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/dwinormalise.log",
    container:
        config["singularity"]["mrtrix"]
    group: "dwiproc"
    shell:
        "dwinormalise individual -nthreads {threads} {input.dwi} {input.mask} {output.dwi} &> {log}"


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
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/dwi2tensor.log",
    group: "dwiproc"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2tensor -nthreads {threads} -mask {input.mask} {input.dwi} {output.dti} &> {log} && "
        "tensor2metric -nthreads {threads} -mask {input.mask} {output.dti} -fa {output.fa} -ad {output.ad} -rd {output.rd} -adc {output.md} >> {log} 2>&1"


# --------------- MRTRIX PREPROC END --------------#


# ------------ MRTRIX TRACTOGRAPHY BEGIN ----------#
rule tckgen:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        fod=rules.mtnormalise.output.wm_fod,
        mask=rules.nii2mif.output.mask,
        cortical_ribbon=rules.fs_xfm_to_native.output.ribbon,
        convex_hull=rules.create_convex_hull.output.convex_hull,
        subcortical_seg=bids_labelmerge(
            space="T1w",
            desc="combined",
            suffix="dseg.nii.gz",
        ),
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
        time=60*18,
    log:
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/tckgen.log",
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tckgen -nthreads {threads} -algorithm iFOD2 -step {params.step} -select {params.sl} -exclude {input.cortical_ribbon} -exclude {input.convex_hull} -include {input.subcortical_seg} -mask {input.mask} -seed_image {input.mask} {input.fod} {output.tck} &> {log}"


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
    threads: 8
    resources:
        mem_mb=32000,
        time=60,
    log:
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/tcksift2.log",
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tcksift2 -nthreads {threads} -out_mu {output.mu} {input.tck} {input.fod} {output.weights} &> {log}"


checkpoint create_roi_mask:
    input:
        subcortical_seg=rules.tckgen.input.subcortical_seg,
        num_labels=rules.get_num_nodes.output.num_labels,
    params:
        base_dir=mrtrix_dir,
        container=config["singularity"]["mrtrix"],
        subj_wildcards=config["subj_wildcards"],
    output:
        out_dir=directory(
            bids_anat_out(
                datatype="roi_masks"
            )
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    group: "tract_masks"
    script:
        "../scripts/mrtpipelines/create_roi_mask.py"


def aggregate_rois(wildcards):
    """Grab all created roi masks"""
    roi_masks = bids_anat_out(
        datatype="roi_masks", 
        desc="{node}", 
        suffix="mask.mif",
        **wildcards,
    )
    # Get node label wildcard
    node = glob_wildcards(roi_masks).node

    # Build node pairs
    node_pairs = np.triu_indices(len(node), k=1)

    return {
        "roi1": expand(
            bids_anat_out(
                datatype="roi_masks",
                desc="{node1}",
                suffix="mask.mif",
            ),
            node1=list(node_pairs[0]+1),
            allow_missing=True,
        ),
        "roi2": expand(
            bids_anat_out(
                datatype="roi_masks",
                desc="{node2}",
                suffix="mask.mif",
            ),
            node2=list(node_pairs[1]+1),
            allow_missing=True,
        )
    }


checkpoint create_exclude_mask:
    input:
        unpack(aggregate_rois),
        rules.create_roi_mask.output.out_dir,
        subcortical_seg=rules.tckgen.input.subcortical_seg,
        num_labels=rules.get_num_nodes.output.num_labels,
    params:
        base_dir=mrtrix_dir,
        mask_dir=bids_anat_out(
            datatype="roi_masks",
        ),
        container=config["singularity"]["mrtrix"],
        subj_wildcards=config["subj_wildcards"],
    output:
        out_dir=directory(
            bids_anat_out(datatype="exclude_mask")
        )
    threads: 4
    resources:
        mem_mb=16000,
        time=60*2,
    group: "tract_masks"
    script:
        "../scripts/mrtpipelines/create_exclude_mask.py"


# TODO (v0.2): ADD OPTION TO OUTPUT TDI MAP


rule tck2connectome:
    """
    Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The effects of SIFT on the reproducibility and biological accuracy of the structural connectome. NeuroImage, 2015, 104, 253-265"

    TODO (v0.2): ADD OPTION TO MULTIPLY BY MU COEFFICIENT
    """
    input:
        weights=rules.tcksift2.output.weights,
        tck=rules.tckgen.output.tck,
        subcortical_seg=rules.tckgen.input.subcortical_seg,
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
        time=60*3,
    log:
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/tck2connectome.log",
    group: "tractography_update"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tck2connectome -nthreads {threads} -zero_diagonal -stat_edge sum -assignment_radial_search {params.radius} -tck_weights_in {input.weights} -out_assignments {output.sl_assignment} -symmetric {input.tck} {input.subcortical_seg} {output.node_weights} &> {log}"


checkpoint connectome2tck:
    input:
        node_weights=rules.tcksift2.output.weights,
        sl_assignment=rules.tck2connectome.output.sl_assignment,
        tck=rules.tckgen.output.tck,
        num_labels=rules.get_num_nodes.output.num_labels,
    params:
        edge_weight_prefix=bids_tractography_out(
            datatype="unfiltered",
            desc="subcortical",
            suffix="tckWeights",
        ),
        edge_tck_prefix=bids_tractography_out(
            datatype="unfiltered",
            desc="subcortical",
            suffix="from",
        ),
    output:
        output_dir=directory(
            bids_tractography_out(
                datatype="unfiltered",
            )
        )
    threads: 32
    resources:
        mem_mb=128000,
        time=60,
    group: "tractography_update"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mkdir -p {output.output_dir} && "
        "num_labels=$(cat {input.num_labels}) && "
        "connectome2tck -nthreads {threads} -nodes `seq -s, 1 $((num_labels))` -exclusive -files per_edge -tck_weights_in {input.node_weights} -prefix_tck_weights_out {params.edge_weight_prefix} {input.tck} {input.sl_assignment} {params.edge_tck_prefix}"


def aggregate_tck_files(wildcards):
    """Grab all files associated with unfiltered tck"""
    # Build list of files
    unfiltered_weights = expand(
        bids_tractography_out(
            datatype="unfiltered",
            desc="subcortical",
            suffix="tckWeights{suffix}.csv",
        ),
        **wildcards,
        allow_missing=True,
    )[0]

    suffix = glob_wildcards(unfiltered_weights).suffix

    unfiltered_tck = expand(
            bids_tractography_out(
                datatype="unfiltered",
                desc="subcortical",
                suffix="from{suffix}.tck",
            ),
            suffix=suffix,
            allow_missing=True,
        )

    return {
        "weights": expand(
            bids_tractography_out(
                datatype="unfiltered",
                desc="subcortical",
                suffix="tckWeights{suffix}.csv",
            ),
            suffix=suffix,
            allow_missing=True
        ),
        "tck": unfiltered_tck, 
    }


def aggregate_exclude_masks(wildcards):
    """Grab all exclude masks"""
    # Build list of files
    exclude_mask = expand(
        bids_anat_out(
            datatype="exclude_mask",
            desc="{desc}",
            suffix="mask.mif",
        ),
        **wildcards,
        allow_missing=True,
    )

    # Get desc wildcard
    desc = glob_wildcards(exclude_mask[0]).desc

    return expand(
        exclude_mask,
        desc=desc,
    )


def get_desc(wildcards):
    """Build param list"""
    # Get desc wildcards
    exclude_mask = expand(
        bids_anat_out(
            datatype="exclude_mask",
            desc="{desc}",
            suffix="mask.mif",
        ),
        **wildcards,
        allow_missing=True,
    )[0]
    desc = glob_wildcards(exclude_mask).desc

    return desc


def get_filtered_tck(wildcards):
    desc = get_desc(wildcards)
    # Create params lists 
    return expand(
        bids_tractography_out(
            desc="{desc}",
            suffix="tractography.tck",
        ),
        **wildcards,
        desc=desc,
        allow_missing=True,
    )


def get_filtered_weights(wildcards):
    desc = get_desc(wildcards)
    return expand(
        bids_tractography_out(
            desc="{desc}",
            suffix="weights.csv",
        ),
        **wildcards,
        desc=desc,
        allow_missing=True,
    )

rule filter_combine_tck:
    input:
        unpack(aggregate_tck_files),
        rules.connectome2tck.output.output_dir,
        rules.create_exclude_mask.output.out_dir,
        filter_mask=aggregate_exclude_masks,
    params:
        filtered_tck=get_filtered_tck,
        filtered_weights=get_filtered_weights,
        filtered_tck_exists=bids_tractography_out(
            desc="from*",
            suffix="tractography.tck"
        ),
        filtered_weights_exists=bids_tractography_out(
            desc="from*",
            suffix="weights.csv"
        ),
        exclude_mask_dir=bids_anat_out(
            datatype="exclude_mask",
        ),
        unfiltered_tck_dir=bids_tractography_out(
            datatype="unfiltered",
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
    threads: 32
    resources:
        mem_mb=128000,
        time=60,
    log:
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/combine_filtered.log",
    group:
        "tractography_update"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "parallel --jobs {threads} -k tckedit -exclude {{1}} -tck_weights_in {{2}} -tck_weights_out {{3}} {{4}} {{5}} ::: {input.filter_mask} :::+ {input.weights} :::+ {params.filtered_weights} :::+ {input.tck} :::+ {params.filtered_tck} || true && " # 'true' to overcome smk bash strict 
        "tckedit {params.filtered_tck_exists} {output.combined_tck} &> {log} && "
        "cat {params.filtered_weights_exists} >> {output.combined_weights} && "
        "rm -r {params.exclude_mask_dir} {params.unfiltered_tck_dir} {params.filtered_tck_exists} {params.filtered_weights_exists}"


rule filtered_tck2connectome:
    input:
        weights=rules.filter_combine_tck.output.combined_weights,
        tck=rules.filter_combine_tck.output.combined_tck,
        subcortical_seg=rules.tckgen.input.subcortical_seg,
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
        time=60*3,
    log:
        f"{config['output_dir']}/logs/mrtrix/sub-{{subject}}/filtered_tck2connectome.log",
    group:
        "tractography_update"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tck2connectome -nthreads {threads} -zero_diagonal -stat_edge sum -assignment_radial_search {params.radius} -tck_weights_in {input.weights} -out_assignments {output.sl_assignment} -symmetric {input.tck} {input.subcortical_seg} {output.node_weights} &> {log}"


# ------------ MRTRIX TRACTOGRAPHY END ----------#

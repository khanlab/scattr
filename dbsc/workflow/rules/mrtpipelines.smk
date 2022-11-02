from os.path import join


include: "zona_bb_subcortex.smk"
include: "freesurfer.smk"


# Directories
mrtrix_dir = join(config["derivatives"], "mrtrix")

# Paramaters
responsemean_flag = config.get("responsemean_dir", None)
shells = config.get("shells", "")
lmaxes = config.get("lmax", "")

# Mrtrix3 citation (additional citations are included per rule as necessary):
# Tournier, J.-D.; Smith, R. E.; Raffelt, D.; Tabbara, R.; Dhollander, T.; Pietsch, M.; Christiaens, D.; Jeurissen, B.; Yeh, C.-H. & Connelly, A. MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. NeuroImage, 2019, 202, 116137


# ------------ MRTRIX PREPROC BEGIN ----------#
rule nii2mif:
    input:
        dwi=inputs["dwi"].input_path,
        bval=lambda wildcards: re.sub(".nii.gz", ".bval", inputs["dwi"].input_path),
        bvec=lambda wildcards: re.sub(".nii.gz", ".bvec", inputs["dwi"].input_path),
        mask=inputs["mask"].input_path,
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
        "mrconvert -nthreads {threads} -fslgrad {input.bvecs} {input.bvals} {input.dwi} {output.dwi} && "
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
        shells=shells,
        lmax=lmaxes,
    output:
        wm_rf=bids(
            root=mrtrix_dir,
            datatype="response",
            desc="wm",
            suffix="response.txt",
            **config["subj_wildcards"],
        ),
        gm_rf=bids(
            root=mrtrix_dir,
            datatype="response",
            desc="gm",
            suffix="response.txt",
            **config["subj_wildcards"],
        ),
        csf_rf=bids(
            root=mrtrix_dir,
            datatype="response",
            desc="csf",
            suffix="response.txt",
            **config["subj_wildcards"],
        ),
    threads: workflow.cores
    group:
        "subject_1"
    container:
        config["singularity"]["mrtrix"]
    shell:
        'if [ {params.shell} == "" ] && [ {params.lmax} == ""]; then'
        "  dwi2response dhollander {input.dwi} {output.wm_rf} {output.gm_rf} {output.csf_rf} -nthreads {threads} -mask {input.mask}"
        "else "
        "  dwi2response dhollander {input.dwi} {output.wm_rf} {output.gm_rf} {output.csf_rf} -nthreads {threads} -mask {input.mask} -shells {params.shells} -lmax {params.lmax} "
        "fi"


rule responsemean:
    """Compute average response function"""
    input:
        wm_rf=expand(rules.dwi2response.output.wm_rf, subject=config["subjects"]),
        gm_rf=expand(rules.dwi2response.output.gm_rf, subject=config["subjects"]),
        csf_rf=expand(rules.dwi2response.output.csf_rf, subject=config["subjects"]),
    output:
        wm_avg_rf=bids(
            root=join(mrtrix_dir, "avg"),
            datatype="response",
            desc="wm",
            suffix="response.txt",
        ),
        gm_avg_rf=bids(
            root=join(mrtrix_dir, "avg"),
            datatype="response",
            desc="gm",
            suffix="response.txt",
        ),
        csf_avg_rf=bids(
            root=join(mrtrix_dir, "avg"),
            datatype="response",
            desc="csf",
            suffix="response.txt",
        ),
    threads: workflow.cores
    group:
        "group"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "responsemean {input.wm_rf} {output.wm_avg_rg} -nthreads {threads} &&"
        "responsemean {input.gm_rf} {output.gm_avg_rg} -nthreads {threads} &&"
        "responsemean {input.csf_rf} {output.csf_avg_rg} -nthreads {threads}"


rule dwi2fod:
    """Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data. NeuroImage, 2014, 103, 411-426"""
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
        wm_rf=join(config["responsemean_dir"], "desc-wm_response.txt")
        if responsemean_flag
        else rules.responsemean.output.wm_avg_rf,
        gm_rf=join(config["responsemean_dir"], "desc-gm_response.txt")
        if responsemean_flag
        else rules.responsemean.output.gm_avg_rf,
        csf_rf=join(config["responsemean_dir"], "desc-csf_response.txt")
        if responsemean_flag
        else rules.responsemean.output.csf_avg_rf,
    params:
        shell=shells,
    output:
        wm_fod=bids(
            root=mrtrix_dir,
            datatype="response",
            model="csd",
            desc="wm",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
        gm_fod=bids(
            root=mrtrix_dir,
            datatype="response",
            model="csd",
            desc="gm",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
        csf_fod=bids(
            root=mrtrix_dir,
            datatype="response",
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
        'if [ {params.shell} == "" ]; then'
        "  dwi2fod -nthreads {threads} -mask {input.mask} msmt_csd {input.dwi} {input.avg_sfwm} {output.wm_fod} {input.avg_gm} {output.gm_fod} {input.avg_csf} {output.avg_csf} "
        "else "
        "  dwi2fod -nthreads {threads} -shell {params.shell} -mask {input.mask} msmt_csd {input.dwi} {input.avg_sfwm} {output.wm_fod} {input.avg_gm} {output.gm_fod} {input.avg_csf} {output.csf} "
        "fi"


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
        wm_fod=bids(
            root=mrtrix_dir,
            datatype="response",
            model="csd",
            desc="normalized",
            suffix="wm_fod.mif",
            **config["subj_wildcards"],
        ),
        gm_fod=bids(
            root=mrtrix_dir,
            datatype="response",
            model="csd",
            desc="normalized",
            suffix="gm_fod.mif",
            **config["subj_wildcards"],
        ),
        csf_fod=bids(
            root=mrtrix_dir,
            datatype="response",
            model="csd",
            desc="normalized",
            suffix="csf_fod.mif",
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
    shell:
        "dwinormalise individual -nthreads {threads} {input.dwi} {input.mask} {output.dwi}"
    group:
        "subject_1"


rule dwi2tensor:
    input:
        dwi=rules.dwinoramlise.dwi,
        mask=rules.nii2mif.output.mask,
    output:
        dti=bids(
            root=mrtrix_dir,
            datatype="dti",
            suffix="dti.mif",
            **config["subj_wildcards"],
        ),
        fa=bids(
            root=mrtrix_dir,
            datatype="dti",
            model="dti",
            suffix="fa.mif",
            **config["subj_wildcards"],
        ),
        ad=bids(
            root=mrtix_dir,
            datatype="dti",
            model="dti",
            suffix="ad.mif",
            **config["subj_wildcards"],
        ),
        rd=bids(
            root=mrtrix_dir,
            datatype="dti",
            model="dti",
            suffix="rd.mif",
            **config["subj_wildcards"],
        ),
        md=bids(
            root=mrtrix_dir,
            datatype="dti",
            model="dti",
            suffix="fa.mif",
            **config["subj_wildcards"],
        ),
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
# TODO (v0.1): CHECK TO MAKE SURE RULES ARE IMPORTED CORRECTLY FROM OTHER SMK FILES
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
        tck=bids(
            root=mrtrix_dir,
            datatype="tractography",
            desc="iFOD2",
            suffix="tractography.tck",
            **config["subj_wildcards"],
        ),
    threads: workflow.cores
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tckgen -nthreads {threads} -algorithm iFOD2 -step {params.step} -select {params.sl_count} -exclude {input.cortical_ribbon} -exclude {input.convex_hull} -include {input.subcortical_seg} -mask {input.mask} -seed_image {input.mask} {input.fod} {output.tck}"


rule tcksift2:
    # Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The effects of SIFT on the reproducibility and biological accuracy of the structural connectome. NeuroImage, 2015, 104, 253-265
    input:
        tck=rules.tckgen.output.tck,
        fod=rules.mtnormalise.output.wm_fod,
    output:
        weights=bids(
            root=mrtrix_dir,
            datatype="tractography",
            desc="iFOD2",
            suffix="tckweights.txt",
            **config["subj_wildcards"],
        ),
        mu=bids(
            root=mrtrix_dir,
            datatype="tractography",
            desc="iFOD2",
            suffix="mucoefficient.txt",
            **config["subj_wildcards"],
        ),
    threads: workflow.cores
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
        sl_assignment=bids(
            root=mrtix_dir,
            datatype="tractography",
            desc="subcortical",
            suffix="nodeassignment.txt",
            **config["subj_wildcards"],
        ),
        node_weights=bids(
            root=mrtrix_dir,
            datatype="tractography",
            desc="subcortical",
            suffix="nodeweights.csv" ** config["subj_wildcards"],
        ),
    threads: workflow.cores
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
    output:
        edge_weight=temp(
            bids(
                root=mrtrix_dir,
                datatype="tractography",
                desc="subcortical",
                suffix="tckweights",
                **config["subj_wildcards"],
            )
        ),
        edge_tck=temp(
            bids(
                root=mrtix_dir,
                datatype="tractography",
                desc="subcortical",
                suffix="from_",
                **config["subj_wildcards"],
            )
        ),
    threads: workflow.cores
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "for i in `seq 2 72`; do "
        "  nodes=$nodes,$i "
        "done "
        "connectome2tck -nthreads {threads} -nodes $nodes -exclusive -filters_per_edge -tck_weights_in {input.node_weights} -prefix_tck_weights_out {output.edge_weight} {input.tck} {input.sl_assignment} {output.edge_tck} "


# NOTE: Use labelmerge split segs here?
rule create_roi_mask:
    # EXPAND OVER NODE1 IN RULE ALL
    input:
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
    output:
        roi_mask=temp(
            bids(
                root=mrtrix_dir,
                datatype="anat",
                desc="{node1}",
                suffix="mask.mif",
                **config["subj_wildcards"]
            )
        ),
    threads: workflow.cores
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrcalc -nthreads {threads} {input.subcortical_seg} {wildcards.node1} -eq {output.roi_mask}"


rule filter_tck:
    # Node1 should be same as prev rule, need to iterate over node2
    input:
        roi1=bids(
            root=mrtrix_dir,
            datatype="anat",
            desc="{node1}",
            suffix="mask.mif",
            **config["subj_wildcards"]
        ),
        roi2=bids(
            root=mrtrix_dir,
            datatype="anat",
            desc="{node2}",
            suffix="mask.mif",
            **config["subj_wildcards"]
        ),
        # ZI rois (do these need to be separately defined in output?)
        lZI=bids(
            root=mrtrix_dir,
            datatype="anat",
            desc="21",
            suffix="mask.mif",
            **config["subj_wildcards"]
        ),
        rZI=bids(
            root=mrtrix_dir,
            datatype="anat",
            desc="22",
            suffix="mask.mif",
            **config["subj_wildcards"]
        ),
        tck=rules.connectome2tck.output.edge_tck,
        weights=rules.tck2connectome.outputs.node_weights,
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
    output:
        filter_mask=temp(
            bids(
                root=mrtrix_dir,
                datatype="anat",
                desc="exclude",
                suffix="mask.mif",
                **config["subj_wildcards"]
            )
        ),
        filtered_tck=temp(
            bids(
                root=mrtrix_dir,
                datatype="tractography",
                desc="from_{node1}-{node2}",
                suffix="tractography.tck",
                **config["subj_wildcards"]
            )
        ),
        filtered_weights=temp(
            bids(
                root=mrtrix_dir,
                datatype="tractography",
                desc="from_{node1}-{node2}",
                suffix="weights.csv",
                **config["subj_wildcards"]
            )
        ),
    threads: workflow.cores
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrcalc -nthreads {threads} {input.subcortical_seg} 0 -neq {input.roi1} -sub {input.roi2} -sub {input.lZI} -sub {input.rZI} -sub {output.filter_mask} "
        "tckedit -nthreads {params.therads} -exclude {output.filter_mask} -tck_weights_in {input.weights} -tck_weights_out {output.filtered_weights} {input.tck} {output.filtered_tck} "


rule combine_filtered:
    input:
        tck=expand(
            rules.filter_tck.output.filtered_tck,
            zip,
            **config["input_zip_lists"]["dwi"]
        ),
        weights=expand(
            rules.filter_tck.output.filtered_weights,
            zip,
            **config["input_zip_lists"]["dwi"]
        ),
    output:
        combined_tck=bids(
            root=mrtrix_dir,
            datatype="tractography",
            desc="subcortical",
            suffix="tractography.tck",
            **config["subj_wildcards"]
        ),
        combined_weights=bids(
            root=mrtrix_dir,
            datatype="tractography",
            desc="subcortical",
            suffix="tckweights.txt",
            **config["subj_wildcards"]
        ),
    threads: workflow.cores
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tckedit {input.tck} {output.combined_tck} "
        "cat {input.weights} >> {output.combined_weights} "


rule filtered_tck2connectome:
    input:
        weights=rules.combine_filtered.output.combined_weights,
        tck=rules.combine_filtered.output.combined_tck,
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
    params:
        radius=config["radial_search"],
    output:
        sl_assignment=bids(
            root=mrtrix_dir,
            datatype="tractography",
            desc="subcortical",
            suffix="nodeassignment.txt",
            **config["subj_wildcards"],
        ),
        node_weights=bids(
            root=mrtrix_dir,
            datatype="tractography",
            desc="subcortical",
            suffix="nodeweights.csv" ** config["subj_wildcards"],
        ),
    threads: workflow.core
    group:
        "subject_2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tck2connectome -nthreads {threads} -zero_diagonal -stat_edge sum -assignment_radial_search {params.radius} -tck_weights_in {input.weights} -out_assignments {output.sl_assignment} -symmetric {input.tck} {input.subcortical_seg} {output.node_weights} -force"

rule nii2mif:
    input:
        dwi=inputs_dwi["dwi"].path,
        bval=re.sub(".nii.gz", ".bval", inputs_dwi["dwi"].path),
        bvec=re.sub(".nii.gz", ".bvec", inputs_dwi["dwi"].path),
        mask=inputs_dwi["mask"].path,
    output:
        dwi=bids(
            root=mrtrix_dir,
            datatype="dwi",
            suffix="dwi.mif",
            **inputs_dwi.subj_wildcards
        ),
        mask=bids(
            root=mrtrix_dir,
            datatype="dwi",
            desc="brain",
            suffix="mask.mif",
            **inputs_dwi.subj_wildcards
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
    log:
        bids_log(suffix="nii2mif.log", **inputs_dwi.subj_wildcards),
    group:
        "dwiproc"
    container:
        config["singularity"]["scattr"]
    conda:
        "../../envs/mrtrix3.yaml"
    shell:
        """
        mrconvert -nthreads {threads} -fslgrad {input.bvec} {input.bval} \\
        {input.dwi} {output.dwi} &> {log} 

        mrconvert -nthreads {threads} {input.mask} {output.mask} >> {log} 2>&1
        """


rule dwi2response:
    """
    Estimate response functions using Dhollander algorithm

    Dhollander, T.; Mito, R.; Raffelt, D. & Connelly, A. 
    Improved white matter response function estimation for 3-tissue 
    constrained spherical deconvolution. Proc Intl Soc Mag Reson Med, 2019, 555
    """
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
    params:
        shells=f"-shells {','.join(shells)}" if shells else "",
        lmax=f"-lmax {','.join(lmax)}" if lmax else "",
        bzero_thresh=config.get("bzero_thresh"),
        mrtrix_conf=temp("~/.mrtrix.conf"),
    output:
        wm_rf=bids_response_out(
            desc="wm",
            **inputs_dwi.subj_wildcards,
        ),
        gm_rf=bids_response_out(
            desc="gm",
            **inputs_dwi.subj_wildcards,
        ),
        csf_rf=bids_response_out(
            desc="csf",
            **inputs_dwi.subj_wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        bids_log(suffix="dwi2response.log"),
    group:
        "dwiproc"
    container:
        config["singularity"]["scattr"]
    conda:
        "../../envs/mrtrix3.yaml"
    shell:
        """
        echo 'BZeroThreshold: {params.bzero_thresh}' > {params.mrtrix_conf}

        dwi2response dhollander {input.dwi} {output.wm_rf} \\
        {output.gm_rf} {output.csf_rf} \\
        -nthreads {threads} -mask {input.mask} \\
        {params.shells} {params.lmax} &> {log}
        """


def get_subject_rf(wildcards):
    """Get appropriate subject response path"""
    if not inputs_dwi.sessions:
        return expand(
            bids_response_out(
                subject="{subject}",
                desc="{tissue}",
            ),
            subject=inputs_dwi.subjects,
            allow_missing=True,
        )
    else:
        return expand(
            bids_response_out(
                subject="{subject}",
                session="{session}",
                desc="{tissue}",
            ),
            subject=inputs_dwi.subjects,
            session=(
                config.get("responsemean_ses")
                if config.get("responsemean_ses")
                else inputs_dwi.sessions
            ),
            allow_missing=True,
        )


rule responsemean:
    """Compute average response function"""
    input:
        subject_rf=get_subject_rf,
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
        f"{log_dir}/{{tissue}}_responsemean.log",
    group:
        "dwiproc_group"
    container:
        config["singularity"]["scattr"]
    conda:
        "../../envs/mrtrix3.yaml"
    shell:
        """
        responsemean {input.subject_rf} {output.avg_rf} \\
        -nthreads {threads} &> {log}
        """


rule dwi2fod:
    """
    Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. 
    Multi-tissue constrained spherical deconvolution for improved analysis 
    of multi-shell diffusion MRI data. NeuroImage, 2014, 103, 411-426
    """
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
        shells=f"-shells {','.join(shells)}" if shells else "",
    output:
        wm_fod=bids_response_out(
            model="csd",
            desc="wm",
            suffix="fod.mif",
            **inputs_dwi.subj_wildcards,
        ),
        gm_fod=bids_response_out(
            model="csd",
            desc="gm",
            suffix="fod.mif",
            **inputs_dwi.subj_wildcards,
        ),
        csf_fod=bids_response_out(
            model="csd",
            desc="csf",
            suffix="fod.mif",
            **inputs_dwi.subj_wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        bids_log(suffix="dwi2fod.log"),
    group:
        "diffmodel"
    container:
        config["singularity"]["scattr"]
    conda:
        "../../envs/mrtrix3.yaml"
    shell:
        """
        dwi2fod -nthreads {threads} {params.shells} -mask {input.mask} \\
        msmt_csd {input.dwi} {input.wm_rf} {output.wm_fod} \\
        {input.gm_rf} {output.gm_fod} \\
        {input.csf_rf} {output.csf_fod} &> {log}
        """


rule mtnormalise:
    """
    Raffelt, D.; Dhollander, T.; Tournier, J.-D.; Tabbara, R.; Smith, R. E.; 
    Pierre, E. & Connelly, A. Bias Field Correction and Intensity 
    Normalisation for Quantitative Analysis of Apparent Fibre Density. 
    In Proc. ISMRM, 2017, 26, 3541

    Dhollander, T.; Tabbara, R.; Rosnarho-Tornstrand, J.; Tournier, J.-D.; 
    Raffelt, D. & Connelly, A. Multi-tissue log-domain intensity and 
    inhomogeneity normalisation for quantitative apparent fibre density. 
    In Proc. ISMRM, 2021, 29, 2472
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
            **inputs_dwi.subj_wildcards,
        ),
        gm_fod=bids_response_out(
            model="csd",
            desc="gm",
            suffix="fodNormalized.mif",
            **inputs_dwi.subj_wildcards,
        ),
        csf_fod=bids_response_out(
            model="csd",
            desc="csf",
            suffix="fodNormalized.mif",
            **inputs_dwi.subj_wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        bids_log(suffix="mtnormalise.log"),
    group:
        "diffmodel"
    container:
        config["singularity"]["scattr"]
    conda:
        "../../envs/mrtrix3.yaml"
    shell:
        """
        mtnormalise -nthreads {threads} -mask {input.mask} \\
        {input.wm_fod} {output.wm_fod} \\
        {input.gm_fod} {output.gm_fod} \\
        {input.csf_fod} {output.csf_fod} &> {log}
        """


rule dwinormalise:
    """DWI normalisation (for DTI)"""
    input:
        dwi=rules.nii2mif.output.dwi,
        mask=rules.nii2mif.output.mask,
    params:
        bzero_thresh=config.get("bzero_thresh"),
        mrtrix_conf=temp("~/.mrtrix.conf"),
    output:
        dwi=bids(
            root=mrtrix_dir,
            datatype="dwi",
            desc="normalized",
            suffix="dwi.mif",
            **inputs_dwi.subj_wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    log:
        bids_log(suffix="dwinormalise.log"),
    container:
        config["singularity"]["scattr"]
    conda:
        "../../envs/mrtrix3.yaml"
    group:
        "dwiproc"
    shell:
        """
        echo 'BZeroThreshold: {params.bzero_thresh}' > {params.mrtrix_conf} 

        dwinormalise individual -nthreads {threads} \\
        {input.dwi} {input.mask} {output.dwi} &> {log}
        """

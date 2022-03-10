from os.path import join

rule nii_to_mif:
    input:
        dwi=bids(
            root=join(config["bids_dir"], 'derivatives/prepdwi')
            datatype='dwi',
            suffix='dwi.nii.gz',
            space='T1w',
            **config['subj_wildcards'],
        )
        bval=bids(
            root=join(config["bids_dir"], 'derivatives/prepdwi')
            datatype='dwi',
            suffix='dwi.bval',
            space='T1w',
            **config['subj_wildcards'],
        )
        bvec=bids(
            root=join(config["bids_dir"], 'derivatives/prepdwi')
            datatype='dwi',
            suffix='dwi.bvec',
            space='T1w',
            **config['subj_wildcards'],
        )
        mask=bids(
            root=join(config["bids_dir"], 'derivatives/prepdwi')
            datatype='dwi',
            suffix='brainmask.nii.gz',
            space='T1w',
            **config['subj_wildcards'],
        )
    params:
        threads=workflow.cores,
    output:
        dwi=bids(
            root=join(config["output_dir"], 'mrtpipelines'),
            datatype='dwi',
            suffix='dwi.mif',
            space='T1w',
            **config['subj_wildcards']
        )
        mask=bids(
            root=join(config['output_dir'], 'mrtpipelines'),
            datatype='dwi',
            suffix='brainmask.mif',
            space='T1w',
            **config['subj_wildcards']
        )
    container:
        config['singularity']['mrtpipelines']
    shell: 
        'mrconvert -nthreads {threads} -fslgrad {input.bvecs} {input.bvals} {input.dwi} {output.dwi} '
        'mrconvert -nthreads {threads} {input.mask} {output.mask}'


rule estimate_response:
    """Estimate response functions using Dhollander algorithm """
    input:
        dwi=rule.nii_to_mif.output.dwi,
        mask=rule.nii_to_mif.output.mask
    params:
        threads=workflow.cores,
    output:
        sfwm=bids(
            root=join(config['output_dir'], 'mrtpipelines'),
            datatype='dwi',
            desc='sfwm',
            suffix='response.txt',
            **config['subj_wildcards'],
        )
        gm=bids(
            root=join(config['output_dir'], 'mrtpipelines'),
            datatype='dwi',
            desc='gm'
            suffix='response.txt',
            **config['subj_wildcards'],
        )
        csf=bids(
            root=join(config['output_dir'], 'mrtpipelines'),
            datatype='dwi',
            desc='csf'
            suffix='response.txt',
            **config['subj_wildcards'],
        )
    container:
        config['singularity']['mrtpipelines']
    shell:
        'dwi2response dhollander {input.dwi} {output.sfwm} {output.gm} {output.csf} -nthreads {threads}'


rule avg_response:
    """Compute average repsonse function"""
    # INPUT NEEDS TO USE ALL SUBJECTS (EXPAND?)
    input:
        sfwm=bids(
                root=join(config['output_dir'], 'mrtpipelines'),
                datatype='dwi',
                desc='sfwm',
                suffix='response.txt',
                **config['subj_wildcards']
            )
        gm=bids(
                root=join(config['output_dir'], 'mrtpipelines'),
                datatype='dwi',
                desc='gm',
                suffix='response.txt',
                **config['subj_wildcards']
            )
        csf=bids(
                root=join(config['output_dir'], 'mrtpipelines'),
                datatype='dwi',
                desc='csf',
                suffix='response.txt',
                **config['subj_wildcards']
            )
    output:
        avg_sfwm=bids(
            root=join(config['output_dir'], 'mrtpipelines/avg'),
            datatype='dwi',
            desc='sfwm',
            suffix='response.txt',
        )
        avg_gm=bids(
            root=join(config['output_dir'], 'mrtpipelines/avg'),
            datatype='dwi',
            desc='gm',
            suffix='response.txt',
        )
        avg_csf=bids(
            root=join(config['output_dir'], 'mrtpipelines/avg'),
            datatype='dwi',
            desc='csf',
            suffix='response.txt',
        )
    container:
        config['singularity']['mrtpipelines']
    shell:
        'average_response {input.sfwm} {output.avg_sfwm} {input.gm} {output.avg_gm} {input.csf} {output.csf}'


if config['shells'] and config['lmax']:
    rule compute_fod:
        input:
            avg_sfwm=rule.avg_response.output.avg_sfwm,
            avg_gm=rule.avg_response.output.avg_gm,
            avg_csf=rule.avg_response.output.avg_csf,
            dwi=rule.nii_to_mif.output.dwi,
            mask=rule.nii_to_mif.output.mask,
        params:
            shell=config['shells'],
            lmax=config['lmax'],
            threads=workflow.threads,
        output:
            wm_fod=bids(
                root=join(config['output_dir'], 'mrtpipelines'),
                datatype='fod',
                model='csd',
                desc='wm',
                suffix='fod.mif'
            ),
            gm_fod=bids(
                root=join(config['output_dir'], 'mrtpipelines'),
                datatype='fod',
                model='csd',
                desc='gm',
                suffix='fod.mif'
            ),
            csf_fod=bids(
                root=join(config['output_dir'], 'mrtpipelines'),
                datatype='fod',
                model='csd',
                desc='csf',
                suffix='fod.mif'
            ),
        container:
            config['singularity']['mrtpipelines']
        shell:
            'dwi2fod -nthreads {params.threads} -shell {params.shell} -lmax {params.lmax} -mask {input.mask} msmt_csd {input.dwi} {input.avg_sfwm} {output.wm_fod} {input.avg_gm} {output.gm_fod} {input.avg_csf} {output.csf}'
else:
    rule compute_fod:
        input:
            avg_sfwm=rule.avg_response.output.avg_sfwm,
            avg_gm=rule.avg_response.output.avg_gm,
            avg_csf=rule.avg_response.output.avg_csf,
            dwi=rule.nii_to_mif.output.dwi,
            mask=rule.nii_to_mif.output.mask,
        params:
            threads=workflow.threads,
        output:
            wm_fod=bids(
                root=join(config['output_dir'], 'mrtpipelines'),
                datatype='fod',
                model='csd',
                desc='wm',
                suffix='fod.mif'
            ),
            gm_fod=bids(
                root=join(config['output_dir'], 'mrtpipelines'),
                datatype='fod',
                model='csd',
                desc='gm',
                suffix='fod.mif'
            ),
            csf_fod=bids(
                root=join(config['output_dir'], 'mrtpipelines'),
                datatype='fod',
                model='csd',
                desc='csf',
                suffix='fod.mif'
            ),
        container:
            config['singularity']['mrtpipelines']
        shell:
            'dwi2fod -nthreads {params.threads} -mask {input.mask} msmt_csd {input.dwi} {input.avg_sfwm} {output.wm_fod} {input.avg_gm} {output.gm_fod} {input.avg_csf} {output.csf}'


rule normalise_fod:
    input:
        wm_fod=rule.compute_fod.output.wm_fod,
        gm_fod=rule.compute_fod.output.gm_fod,
        csf_fod=rule.compute_fod.output.csf_fod,
        mask=rule.nii_to_mif.output.mask,
    params:
        threads=workflow.threads,
    output:
        wm_fod=bids(
            root=join(config['output_dir'], 'mrtpipelines'),
            datatype='fod',
            model='csd',
            desc='wm',
            suffix='fodnorm.mif'
        ),
        gm_fod=bids(
            root=join(config['output_dir'], 'mrtpipelines'),
            datatype='fod',
            model='csd',
            desc='gm',
            suffix='fodnorm.mif'
        ),
        csf_fod=bids(
            root=join(config['output_dir'], 'mrtpipelines'),
            datatype='fod',
            model='csd',
            desc='csf',
            suffix='fodnorm.mif'
        ),
    container:
        config['singularity']['mrtpipelines']
    shell:
        'mtnormalise -nthreads {params.threads} -mask {input.mask} {input.wm_fod} {output.wm_fod} {input.gm_fod} {output.gm_fod} {input.csf_fod} {output.csf_fod}'

## ADD TENSOR PROCESSING
## ADD TRACTOGRAPHY
## ADD SUBCORTICAL STUFF
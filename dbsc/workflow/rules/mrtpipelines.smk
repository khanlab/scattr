from os.path import join

include: 'zona_bb_subcortex.smk'
include: 'freesurfer.smk'

# DWI Processing
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
        dwi=rules.nii_to_mif.output.dwi,
        mask=rules.nii_to_mif.output.mask
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
    """Compute average response function"""
    input:
        sfwm=expand(rules.estimate_response.output.sfwm, zip, **config['input_zip_lists']['T1w']),
        gm=expand(rules.estimate_response.output.gm, zip, **config['input_zip_lists']['T1w']),
        csf=expand(rules.estimate_response.output.csf, zip, **config['input_zip_lists']['T1w'])
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


rule compute_fod:
    input:
        avg_sfwm=rules.avg_response.output.avg_sfwm,
        avg_gm=rules.avg_response.output.avg_gm,
        avg_csf=rules.avg_response.output.avg_csf,
        dwi=rules.nii_to_mif.output.dwi,
        mask=rules.nii_to_mif.output.mask,
    params:
        shell=config.get('shells', ''),
        lmax=config.get('lmax', ''),
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
        'if [ {params.shell} == "" ] && [ {params.lmax} == "" ]; then'
        '  dwi2fod -nthreads {params.threads} -shell {params.shell} -lmax {params.lmax} -mask {input.mask} msmt_csd {input.dwi} {input.avg_sfwm} {output.wm_fod} {input.avg_gm} {output.gm_fod} {input.avg_csf} {output.csf} '
        'else '
        '  dwi2fod -nthreads {params.threads} -mask {input.mask} msmt_csd {input.dwi} {input.avg_sfwm} {output.wm_fod} {input.avg_gm} {output.gm_fod} {input.avg_csf} {output.avg_csf} '
        'fi'


rule normalise_fod:
    input:
        wm_fod=rules.compute_fod.output.wm_fod,
        gm_fod=rules.compute_fod.output.gm_fod,
        csf_fod=rules.compute_fod.output.csf_fod,
        mask=rules.nii_to_mif.output.mask,
    params:
        threads=workflow.threads,
    output:
        wm_fod=bids(
            root=join(config['output_dir'], 'mrtpipelines'),
            datatype='fod',
            model='csd',
            desc='wm',
            suffix='fodnorm.mif',
            **config['subj_wildcards'],
        ),
        gm_fod=bids(
            root=join(config['output_dir'], 'mrtpipelines'),
            datatype='fod',
            model='csd',
            desc='gm',
            suffix='fodnorm.mif',
            **config['subj_wildcards'],
        ),
        csf_fod=bids(
            root=join(config['output_dir'], 'mrtpipelines'),
            datatype='fod',
            model='csd',
            desc='csf',
            suffix='fodnorm.mif',
            **config['subj_wildcards'],
        ),
    container:
        config['singularity']['mrtpipelines']
    shell:
        'mtnormalise -nthreads {params.threads} -mask {input.mask} {input.wm_fod} {output.wm_fod} {input.gm_fod} {output.gm_fod} {input.csf_fod} {output.csf_fod}'


# DTI (Tensor) Processing
rule dwi_normalise:
    input:
        dwi=rules.nii_to_mif.output.dwi,
        mask=rules.nii_to_mif.output.mask,
    params:
        threads=workflow.threads,
    output:
        dwi=bids(
            root=join(config['out_dir'], 'mrtpipelines'),
            datatype='dwi',
            space='T1w',
            desc='norm',
            suffix='dwi.mif',
            **config['subj_wildcards'],
        )
    container:
        config['singularity']['mrtpipelines']
    shell:
        'dwinormalise -nthreads {params.threads} {input.dwi} {input.mask} {output.dwi}'


rule compute_tensor:
    input:
        dwi=rules.dwi_noramlise.dwi,
        mask=rules.nii_to_mif.output.mask,
    params:
        threads=workflow.threads,
    output:
        dti=bids(
            root=join(config['out_dir'], 'mrtpipelines'),
            datatype='dti',
            space='T1w',
            suffix='dti.mif',
            **config['subj_wildcards'],
        ),
        fa=bids(
            root=join(config['out_dir'], 'mrtpipelines'),
            datatype='dti',
            space='T1w',
            model='dti',
            suffix='fa.mif',
            **config['subj_wildcards'],
        ),
        ad=bids(
            root=join(config['out_dir'], 'mrtpipelines'),
            datatype='dti',
            space='T1w',
            model='dti',
            suffix='ad.mif',
            **config['subj_wildcards'],
        ),
        rd=bids(
            root=join(config['out_dir'], 'mrtpipelines'),
            datatype='dti',
            space='T1w',
            model='dti',
            suffix='rd.mif',
            **config['subj_wildcards'],
        ),
        md=bids(
            root=join(config['out_dir'], 'mrtpipelines'),
            datatype='dti',
            space='T1w',
            model='dti',
            suffix='fa.mif',
            **config['subj_wildcards'],
        )
    container:
        config['singularity']['mrtpipelines']
    shell:
        'dwi2tensor -nthreads {params.threads} -mask {input.mask} {input.dwi} {output.dti}'
        'tensor2metric -nthreads {params.threads} -mask {input.mask} {output.dti} -fa {output.fa} -ad {output.ad} -rd {output.rd} -adc {output.md}'


# Tractography processing
rule gen_tractography:
    input:
        fod=rules.normalise_fod.output.wm_fod,
        cortical_ribbon=#ADD CORTICAL RIBBON
        convex_hull=rules.create_convex_hull.output.convex_hull,
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
        mask=rules.nii_to_mif.output.mask,
    params:
        threads=workflow.threads,
        step=config['step'],
        sl=config['sl_count'],
    output:
        tck=bids(
            root=join(config['out_dir'], mrtpipelines),
            datatype='tractography',
            space='T1w',
            desc='iFOD2',
            suffix='tractography.tck',
            **config["subj_wildcards"],
        )
    container:
        config['singularity']['mrtpipelines']
    shell:
        'tckgen -nthreads {params.threads} -algorithm iFOD2 -step {params.step} -select {params.sl_count} -exclude {input.cortical_ribbon} -exclude {input.convex_hull} -include {input.subcortical_seg} -mask {input.mask} -seed_image {input.mask} {input.fod} {output.tck}'

rule weight_tractography:
    input:
        tck=rules.gen_tractography.output.tck,
        fod=rules.normalise_fod.output.wm_fod,
    params:
        threads=workflow.threads,
    output:
        weights=bids(
            root=join(config['out_dir'], 'mrtpipelines'),
            datatype='tractography',
            space='T1w',
            desc='iFOD2',
            suffix='tckweights.txt',
            **config['subj_wildcards'],
        ),
        mu=bids(
            root=join(config['out_dir'], 'mrtpipelines'),
            datatype='tractography',
            space='T1w',
            desc='iFOD2',
            suffix='mucoefficient.txt',
            **config['subj_wildcards'],
            ]
        ),
    container:
        config['singularity']['mrtpipelines']
    shell:
        'tcksift2 -nthreads {params.threads} -out_mu {output.mu} {input.tck} {input.fod} {output.weights}'

# ADD OPTION TO OUTPUT TDI MAP

# Connectivity map
# ADD OPTION TO MULTIPLY BY MU COEFFICIENT
rule connectome_map:
    input:
        weights=rules.weight_tractography.output.weights,
        tck=rules.gen_tractography.output.tck, 
        subcortical_seg=rules.add_brainstem_new_seg.output.seg,
    params:
        threads=workflow.threads,
        radius=config['radial_search']
    output:
        sl_assignment=bids(
            root=join(config['out_dir'], 'mrtpipelines'),
            datatype='tractography',
            space='T1w',
            desc='subcortical',
            suffix='nodeassignment.txt',
            **config['subj_wildcards'],
        ),
        node_weights=bids(
            root=join(config['out_dir'], 'mrtpipelines'),
            datatype='tractography',
            space='T1w',
            desc='subcortical',
            suffix='nodeweights.csv'
            **config['subj_wildcards'],
        )
    container:
        config['singularity']['mrtpipelines']
    shell:
        'tck2connectome -nthreads {params.threads} -zero_diagonal -stat_edge sum -assignment_radial_search {params.radius} -tck_weights_in {input.weights} -out_assignments {output.sl_assignment} -symmetric {input.tck} {input.subcortical_seg} {output.node_weights} '

# Remove streamlines passing through other GM ROI
if config['exclude_gm']:
    rule extract_tck:
        input:
            node_weights=rules.connectome_map.output.node_weights,
            sl_assignment=rules.connectome_map.output.sl_assignment,
            tck=rules.gen_tractography.output.tck, 
        params:
            threads=workflow.threads,
        output:
            edge_weight=temp(
                bids(
                    root=join(config['out_dir'], 'mrtpipelines'),
                    datatype='tractography',
                    space='T1w',
                    desc='subcortical',
                    suffix='tckweights',
                    **config['subj_wildcards'],
                )
            ),
            edge_tck=temp(
                bids(
                    root=join(config['out_dir'], 'mrtpipelines'),
                    datatype='tractography',
                    space='T1w',
                    desc='subcortical',
                    suffix='from_',
                    **config['subj_wildcards'],
                )
            ),
        container:
            config['singularity']['mrtpipelines'],
        shell:
            'for i in `seq 2 72`; do '
            '  nodes=$nodes,$i '
            'done '
            'connectome2tck -nthreads {params.threads} -nodes $nodes -exclusive -filters_per_edge -tck_weights_in {input.node_weights} -prefix_tck_weights_out {output.edge_weight} {input.tck} {input.sl_assignment} {output.edge_tck} '


    rule create_roi_mask:
        # EXPAND OVER NODE1 IN RULE ALL
        input:
            subcortical_seg=rules.add_brainstem_new_seg.output.seg,
        params:
            threads=workflow.threads,
        output:
            roi_mask=temp(
                bids(
                    root=join(config['out_dir'], 'mrtpipelines'),
                    datatype='anat',
                    space='T1w',
                    desc='{node1}',
                    suffix='mask.mif'
                )
            ),
        container:
            config['singularity']['mrtpipelines'],
        shell:
            'mrcalc -nthreads {params.threads} {input.subcortical_seg} {wildcards.node1} -eq {output.roi_mask}'

    rule filter_tck:
        # Node1 should be same as prev rule, need to iterate over node2
        input:
            roi1=bids(
                root=join(config['out_dir'], 'mrtpipelines'),
                datatype='anat',
                space='T1w',
                desc='{node1}',
                suffix='mask.mif'
            ),
            roi2=bids(
                root=join(config['out_dir'], 'mrtpipelines'),
                datatype='anat',
                space='T1w',
                desc='{node2}',
                suffix='mask.mif'
            ),
            # ZI rois (do these need to be separately defined in output?)
            lZI=bids(
                root=join(config['out_dir'], 'mrtpipelines'),
                datatype='anat',
                space='T1w',
                desc='21',
                suffix='mask.mif'
            ),            
            rZI=bids(
                root=join(config['out_dir'], 'mrtpipelines'),
                datatype='anat',
                space='T1w',
                desc='22',
                suffix='mask.mif'
            ),
            subcortical_seg=rules.add_brainstem_new_seg.output.seg,
            tck=rules.extract_tck.output.edge_tck,
            weights=rules.connectome_map.outputs.node_weights
        params:
            threads=workflow.threads,
        output:
            filter_mask=temp(
                bids(
                    root=join(config['out_dir'], 'mrtpipelines'),
                    datatype='anat',
                    space='T1w',
                    desc='exclude',
                    suffix='mask.mif'
                )
            )
            filtered_tck=temp(
                    bids(
                root=join(config['out_dir'], 'mrtpipelines'),
                datatype='tractography',
                space='T1w',
                desc='from_{node1}-{node2}',
                suffix='tractography.tck'
                )
            )
            filtered_weights=temp(
                    bids(
                    root=join(config['out_dir'], 'mrtpipelines'),
                    datatype='tractography',
                    space='T1w',
                    desc='from_{node1}-{node2}',
                    suffix='weights.csv'
                )
            )
        containers:
            config['singularity']['mrtpipelines'],
        shell: 
            'mrcalc -nthreads {params.threads} {input.subcortical_seg} 0 -neq {input.roi1} -sub {input.roi2} -sub {input.lZI} -sub {input.rZI} -sub {output.filter_mask} '
            'tckedit -nthreads {params.therads} -exclude {output.filter_mask} -tck_weights_in {input.weights} -tck_weights_out {output.filtered_weights} {input.tck} {output.filtered_tck} '


    rule combine_filtered:
        input:
            tck=expand(rules.filter_tck.output.filtered_tck, zip, **config['input_zip_lists']['T1w']),
            weights=expand(rules.filter_tck.output.filtered_weights, zip, **config['input_zip_lists']['T1w']),
        params:
            threads=workflow.threads,
        output:
            combined_tck=bids(
                root=join(config['out_dir'], 'mrtpipelines'),
                datatype='tractography',
                space='T1w',
                desc='subcortical',
                suffix='tractography.tck'
            )
            combined_weights=bids(
                root=join(config['out_dir'], 'mrtpipelines'),
                datatype='tractography',
                space='T1w',
                desc='subcortical',
                suffix='tckweights.txt'
            )
        container:
            config['singularity']['mrtpipelines'],
        shell:
            'tckedit {input.tck} {output.combined_tck} '
            'cat {input.weights} >> {output.combined_weights} '

    rule filtered_connectome_map:
        input:
            weights=rules.combine_filtered.output.combined_weights,
            tck=rules.combine_filtered.output.combined_tck, 
            subcortical_seg=rules.add_brainstem_new_seg.output.seg,
        params:
            threads=workflow.threads,
            radius=config['radial_search']
        output:
            sl_assignment=bids(
                root=join(config['out_dir'], 'mrtpipelines'),
                datatype='tractography',
                space='T1w',
                desc='subcortical',
                suffix='nodeassignment.txt',
                **config['subj_wildcards'],
            ),
            node_weights=bids(
                root=join(config['out_dir'], 'mrtpipelines'),
                datatype='tractography',
                space='T1w',
                desc='subcortical',
                suffix='nodeweights.csv'
                **config['subj_wildcards'],
            )
        container:
            config['singularity']['mrtpipelines']
        shell:
            'tck2connectome -nthreads {params.threads} -zero_diagonal -stat_edge sum -assignment_radial_search {params.radius} -tck_weights_in {input.weights} -out_assignments {output.sl_assignment} -symmetric {input.tck} {input.subcortical_seg} {output.node_weights} -force'

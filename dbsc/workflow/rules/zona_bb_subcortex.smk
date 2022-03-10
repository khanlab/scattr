from os.path import join 

rule xfm_to_native:
    input:
        seg=config['zona_bb_subcortex'][config['Space']]['seg']
        ref=bids(
            root=config["bids_dir"],
            datatype="anat",
            suffix="T1w.nii.gz",
            **config["subj_wildcards"],
        ),
        xfm=## ADD XFM from template to native
    output:
        nii=bids(
            root=join(config['output_dir'], 'zona_bb_subcortex'),
            datatype="anat",
            space='T1w',
            desc='ZonaBB'
            suffix="seg.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config['singularity']['neuroglia-core'],
    shell: 
        'applywarp --rel --interp=nn -i {input.seg} -r {input.ref} -w {input.xfm} -o {output.nii}'


rule binarize:
    input: 
        nii=rule.xfm_to_native.output.nii,
    output:
        bin=bids(
            root=join(config['output_dir'], 'zona_bb_subcortex'),
            datatype="anat",
            space='T1w',
            desc='ZonaBB'
            suffix="seg_bin.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config['singularity']['neuroglia-core'],
    shell:
        'fslmaths {input.nii} -bin {output.bin}'


rule add_brainstem:
    ## NEED TO UPDATE BRAINSTEM INPUT
    input:
        bin=rule.binarize.output.bin,
        aparcaseg=#APARCASEG.NII.GZ
    output:
        bin=bids(
            root=join(config['output_dir'], 'zona_bb_subcortex'),
            datatype="anat",
            space='T1w',
            desc='ZonaBBStem'
            suffix="seg_bin.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config['singularity']['neuroglia-core']
    shell:
        'fslmaths {input.aparcaseg} -thr 16 -uthr 16 -bin -max {input.bin} {output.bin}'


rule xfm_zona_rois:
    input: 
        mask=expand(
            join(
                config['zona_bb_subcortex'][config['Space']]['dir'],
                f'sub-SNSX32Nlin2020Asym_space-{config["Space"]}_hemi-{{hemi}}_desc-{{struct}}_mask.nii.gz'
            )
            hemi=['L', 'R'],
            struct=['fct', 'ft', 'fl', 'hfields']
        )
        ref=rule.xfm_to_native.input.ref,
        xfm=rule.xfm_to_native.input.xfm,
    output:
        mask=expand(
            bids(
                root=join(config['output_dir'], 'zona_bb_subcortex'),
                datatype="anat",
                space="T1w",
                hemi="{hemi}",
                desc="{struct}",
                suffix="mask.nii.gz",
            )
            hemi=['L', 'R'],
            struct=['fct', 'ft', 'fl', 'hfields']
        )
    container:
        config['singularity']['neuroglia-core']
    shell:
        'applywarp --rel --interp=nn -i {input.mask} -r {input.ref} -w {input.xfm} -o {output.mask}'


rule rm_bb_thal:
    """Removes existing thalamus"""
    input:
        seg=rule.xfm_to_native.output.nii, 
    output:
        seg=bids(
            root=join(config['output_dir'], 'zona_bb_subcortex'),
            datatype="anat",
            space='T1w',
            desc='ZonaBBNoThal',
            suffix="seg.nii.gz",
            **config["subj_wildcards"],
    ),
        non_thal=temp(
            bids(
                root=join(config['output_dir'], 'zona_bb_subcortex'),
                datatype="anat",
                space='T1w',
                desc='NonThal',
                suffix="seg.nii.gz",
                **config["subj_wildcards"], 
            )
        )
        rm_seg=temp(
            bids(
                root=join(config['output_dir'], 'zona_bb_subcortex'),
                datatype="anat",
                space='T1w',
                desc='rm',
                suffix="seg.nii.gz",
                **config["subj_wildcards"], 
            )
        )
    container:
        config['singularity']['neuroglia-core']
    shell:
        'fslmaths {input.seg} -sub 2 {output.non_thal} '
        'fslmaths {output.non_thal} -thr 15 {output.non_thal} '
        'fslmaths {input.seg} -thr 15 {output.rm_seg} '
        'fslmaths {input.seg} -sub {output.rm_seg} {output.seg}'
        'fslmaths {output.seg} -add {output.non_thal} {output.seg}'


rule add_new_thal:
    input: 
        aparcaseg=rule.add_brainstem.input.aparcaseg, 
        labels=config['freesurfer']['labels'],
        thal=rule.freesurfer.xfm_to_native.output.nii,
        seg=rule.rm_bb_thal.output.rm_seg,
    output:
        seg=bids(
            root=join(config['output_dir'], 'zona_bb_subcortex'),
            datatype="anat",
            space='T1w',
            desc='ZonaBBFSThal',
            suffix="seg.nii.gz",
            **config["subj_wildcards"],
        ),
        label=temp(
            bids(
                root=join(config['output_dir'], 'zona_bb_subcortex'),
                datatype='anat',
                space='T1w',
                suffix='label.nii.gz',
                **config['subj_wildcards'],
            )
        )
    container:
        config['singularity']['neuroglia-core']
    shell:
        'labelIdx=23 '
        'cp {input.seg} {output.seg} '
        'while IFS=, read -r label hemi nuclei || [ -n "$label" ]; do '
        '  if [ "$label" = "# Label" ]; then '
        '    continue '
        '  fi '
        '  fslmaths {input.thal} -thr $label -uthr $label -bin {output.label} '
        '  fslmaths {output.label} -mul $labelIdx {output.label} '
        '  fslmaths {output.seg} -max {output.label} {output.seg} '
        '  labelIdx=$((labelIdx+1)) '
        ' done < {input.labels}'


rule bin_new_seg:
    input:
        seg=rule.add_new_thal.output.seg
    output:
        seg=bids(
            root=join(config['output_dir'], 'zona_bb_subcortex'),
            datatype="anat",
            space='T1w',
            desc='ZonaBBFSThal',
            suffix="seg_bin.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config['singularity']['neuroglia-core']
    shell:
        'fslmaths {input.seg} -bin {output.seg}'


rule add_brainstem_new_seg:
    input:
        seg=rule.bin_new_seg.output.seg
        aparcaseg=rule.add_brainstem.input.aparcaseg
    output:
        seg=bids(
            root=join(config['output_dir'], 'zona_bb_subcortex'),
            datatype='anat',
            space='T1w',
            desc='ZonaBBFSThalStem',
            suffix='seg_bin.nii.gz',
            **config['subj_wildcards'],
        )
    container:
        config['singularity']['neuroglia-core']
    shell:
        'fslmaths {input.aparcaseg} -thr 16 -uthr 16 -bin -max {input.seg} {output.seg}'
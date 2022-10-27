from os.path import join 

include: "freesurfer.smk"

# ADD SEGMENTATION FOR CORTICAL RIBBON

rule xfm_to_native:
    input:
        seg=config['zona_bb_subcortex'][config['Space']]['seg']
        ref=bids(
            root=config["bids_dir"],
            datatype="anat",
            suffix="T1w.nii.gz",
            **config["subj_wildcards"],
        ),
        xfm=#ADD XFM from template to native#
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
        nii=rules.xfm_to_native.output.nii,
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
    input:
        bin=rules.binarize.output.bin,
        aparcaseg=rules.fs_xfm_to_native.output.aparcaseg,
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
        ref=rules.xfm_to_native.input.ref,
        xfm=rules.xfm_to_native.input.xfm,
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

# Everything in the block below should be replaced by labelmerge
######################################################################
rule rm_bb_thal:
    """Removes existing thalamus"""
    input:
        seg=rules.xfm_to_native.output.nii, 
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
        aparcaseg=rules.add_brainstem.input.aparcaseg, 
        labels=config['freesurfer']['labels'],
        thal=rules.fs_xfm_to_native.output.thal,
        seg=rules.rm_bb_thal.output.rm_seg,
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
        seg=rules.add_new_thal.output.seg
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
        seg=rules.bin_new_seg.output.seg
        aparcaseg=rules.add_brainstem.input.aparcaseg
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
###############################################################################

rule create_convex_hull:
    input:
        bin_seg=rules.add_brainstem_new_seg.output.seg,
    output:
        convex_hull=bids(
            root=join(config['output_dir'], 'zona_bb_subcortex'),
            datatype='anat',
            space='T1w',
            desc='ConvexHull',
            suffix=seg_bin.nii.gz',
            **config['subj_wildcards'],
        )
    script:
        '../resources/zona_bb_subcortex/convexHull_roi.py'
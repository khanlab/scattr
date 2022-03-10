#!/bin/bash

# User Inputs
subj=sub-${1:4}

# Directories
ZONA_DIR=$HOME/scratch/Zona
OUT_DIR=$ZONA_DIR/data/hcpur100_3T/derivatives/zona_bb_subcortex/$subj/anat
localscratch=$SLURM_TMPDIR

# Files / Containers
SINGULARITY_IMG=/home/tkai/singularity/khanlab_neuroglia-core_v1.3.0.img
BBSEG=$ZONA_DIR/misc/zona_bb_subcortex/output_MNI2009b_to_MNI152/ZonaBBSubCorSeg_MNI152.nii.gz
REF_IMG=$ZONA_DIR/data/hcpur100_3T/$subj/anat/${subj}_acq-procHCP_T1w.nii.gz
XFM_WARP=$localscratch/$subj/xfms/standard2acpc_dc.nii.gz
FS_SEG=$ZONA_DIR/data/hcpur100_3T/derivatives/freesurfer/$subj/anat/${subj}_space-T1w_aparc+aseg.nii.gz

# Scripts
mkdir -p $OUT_DIR

## Unzip xfm
mkdir -p $localscratch/$subj/xfms 

unzip $HOME/projects/ctb-akhanf/ext-data/hcp1200/zipfiles/${1:4}_3T_Structural_preproc.zip ${1:4}/MNINonLinear/xfms/standard2acpc_dc.nii.gz -d $localscratch/$subj/xfms
mv $localscratch/$subj/xfms/${1:4}/MNINonLinear/xfms/standard2acpc_dc.nii.gz $localscratch/$subj/xfms/standard2acpc_dc.nii.gz

## Transform subcortical segmentation
singularity exec $SINGULARITY_IMG applywarp --rel --interp=nn -i $BBSEG -r $REF_IMG -w $XFM_WARP -o $OUT_DIR/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz

singularity exec $SINGULARITY_IMG fslmaths $OUT_DIR/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz -bin $OUT_DIR/${subj}_space-T1w_desc-ZonaBBSubCorSeg_bin.nii.gz

#Add FS brainstem
echo "Adding brainstem to binary mask"
singularity exec $SINGULARITY_IMG fslmaths $FS_SEG -thr 16 -uthr 16 -bin -max $OUT_DIR/${subj}_space-T1w_desc-ZonaBBSubCorSeg_bin.nii.gz $OUT_DIR/${subj}_space-T1w_desc-ZonaBBSubCorSeg_bin.nii.gz


## Transform zona ROIs
for hemi in R L; do
  for struct in fl ft fct hfields; do 
    singularity exec $SINGULARITY_IMG applywarp --rel --interp=nn -i $ZONA_DIR/misc/zona_bb_subcortex/output_MNI2009b_to_MNI152/MNI152_sub-SNSX32NLin2020Asym_hemi-${hemi}_desc-${struct}_mask.nii.gz -r $REF_IMG -w $XFM_WARP -o $OUT_DIR/${subj}_space-T1w_hemi-${hemi}_desc-${struct}_mask.nii.gz
  done
done


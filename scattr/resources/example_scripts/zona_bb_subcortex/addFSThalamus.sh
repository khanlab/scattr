#!/bin/bash
# Setup
subj=$1 
njobs=$2
localscratch=$SLURM_TMPDIR
ZONA_DIR=$HOME/scratch/Zona/data/hcpur100_3T
THAL_PARC=$ZONA_DIR/derivatives/freesurfer/$subj/anat/${subj}_space-T1w_desc-fs_thalamus.nii.gz
FS_SEG=$ZONA_DIR/derivatives/freesurfer/$subj/anat/${subj}_space-T1w_aparc+aseg.nii.gz

SINGULARITY_IMG=/home/tkai/singularity/khanlab_neuroglia-core_v1.3.0.img
MRTRIX_IMG=/home/tkai/opt/singularity/mrtrix3-dev.sif
SINGULARITYENV_OMP_NUM_THREADS=$njobs
SINGULARITYENV_MKL_NUM_THREADS=$njobs

# Processing
mkdir -p $localscratch/$subj/anat

## Old thalamus
echo "Removing old thalamus labels and updating..."
cp $ZONA_DIR/derivatives/zona_bb_subcortex/${subj}/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz

### Grab labels following thalamus
singularity exec $SINGULARITY_IMG fslmaths $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz -sub 2 $localscratch/$subj/anat/nonthal_labels.nii.gz

singularity exec $SINGULARITY_IMG fslmaths $localscratch/$subj/anat/nonthal_labels.nii.gz -thr 15 $localscratch/$subj/anat/nonthal_labels.nii.gz

### Remove existing thalamus labels
singularity exec $SINGULARITY_IMG fslmaths $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz -thr 15 $localscratch/$subj/anat/rm_labels.nii.gz

singularity exec $SINGULARITY_IMG fslmaths $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz -sub $localscratch/$subj/anat/rm_labels.nii.gz $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz

### Update labels following thalamus
singularity exec $SINGULARITY_IMG fslmaths $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz -add $localscratch/$subj/anat/nonthal_labels.nii.gz $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz


## New thalamus
echo "Adding thalamus parcellation from Freesurfer 7.x..."
labelIdx=23
while IFS=, read -r label hemi nuclei || [ -n "$label" ]; do
  if [ "$label" = "# Label" ]; then
    continue
  fi
  echo "Currently adding: $hemi $nuclei..."

  # Extract label
  singularity exec $SINGULARITY_IMG fslmaths $THAL_PARC -thr $label -uthr $label -bin $localscratch/$subj/anat/label.nii.gz

  singularity exec $SINGULARITY_IMG fslmaths $localscratch/$subj/anat/label.nii.gz -mul $labelIdx $localscratch/$subj/anat/label.nii.gz

  # Add to parcellation
  singularity exec $SINGULARITY_IMG fslmaths $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz -max $localscratch/$subj/anat/label.nii.gz $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz

  # Increment labelIdx
  labelIdx=$((labelIdx+1))
done < $ZONA_DIR/derivatives/freesurfer/code/fs_labels.csv

echo "Recreating binarization of mask"
singularity exec $SINGULARITY_IMG fslmaths $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz -bin $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg_bin.nii.gz

#Add FS brainstem
echo "Adding brainstem to binary mask"
singularity exec $SINGULARITY_IMG fslmaths $FS_SEG -thr 16 -uthr 16 -bin -max $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg_bin.nii.gz $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg_bin.nii.gz

# Syncing 
echo "Syncing data..."
rsync --update $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz $ZONA_DIR/derivatives/zona_bb_subcortex/${subj}/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg.nii.gz
rsync --update $localscratch/$subj/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg_bin.nii.gz $ZONA_DIR/derivatives/zona_bb_subcortex/${subj}/anat/${subj}_space-T1w_desc-ZonaBBSubCorSeg_bin.nii.gz

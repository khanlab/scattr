## ANTS transforms from MNI2009b to MNI152 (HCP)

```
antsApplyTransforms -d 3 -i <input_space-MNI2009b> -o <output_space-MNI152> -r <reference_space-MNI152> -t ./mni2009b-to-mni152-T1w1Warp.nii.gz [./mni2009b-to-mni152-T1w0GenericAffine.mat,0]
```

`input_space-MNI2009b` - input image in MNI2009bAsym space
`output_space-MNI152` - output image in MNI152NLin6Asym space
`<reference_space-MNI152` - reference image for spacing, voxel size, etc.

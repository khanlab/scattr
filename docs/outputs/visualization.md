# Visualization

## MRview

[MRview](https://mrtrix.readthedocs.io/en/latest/reference/commands/mrview.html)
is a powerful viewer distributed with [MRtrix](https://www.mrtrix.org/) 
that works well with both volume (e.g. `.nii`, `.nii.gz`) and tractography
(e.g. `.tck`) data. It may take a while to load and navigate large 
tractography files.

### Tips and tricks:
* Consider extracting a particular connection passing through two anatomical 
structures to view instead of the full tractogram.  This can be done with the 
tractogram data, and the structures of interest:

```
tckedit -include structure1_roi.nii.gz -include structure2_roi.nii.gz in_tractogram.tck out_tractogram.tck
```

* Extracting a subset number of streamlines from tracts (e.g. 10K):

```
tckedit in_tractogram.tck out_tractogram.tck -number 1000
```

_Note: MRtrix has a number of different options for manipulating tractograms 
which can aid the ability to visualize such tracts. Take a look through the 
[`tckedit` documentation](https://mrtrix.readthedocs.io/en/latest/reference/commands/tckedit.html)
as well as both the [`tck2connectome`](https://mrtrix.readthedocs.io/en/latest/reference/commands/tck2connectome.html) and 
[`connectome2tck` documentation](https://mrtrix.readthedocs.io/en/latest/reference/commands/connectome2tck.html)._

## ITK-SNAP

ITK-SNAP is a lightweight tool that is able to quickly open volumes and is ideal
for viewing segmentation images (`_dseg.nii.gz`). Segmentations can be loaded
as overlays and ITK-SNAP can create a 3D rendering of the contours of each
label. 
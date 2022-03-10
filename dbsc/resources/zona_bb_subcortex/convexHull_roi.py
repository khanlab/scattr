#!/bin/env python
import nibabel as nib 
import numpy as np 
import scipy.spatial 

def create_convex_hull(binary_seg, out_convex_hull):
    # Load input  
    print("Loading data...")   
    seg = nib.load(binary_seg)
    seg_affine = seg.affine
    seg_data = seg.get_fdata()

    # Compute convex hull
    print("Computing convex hull...")
    points = np.transpose(np.where(seg_data))
    hull = scipy.spatial.ConvexHull(points)
    deln = scipy.spatial.Delaunay(points[hull.vertices])
    idx = np.stack(np.indices(seg_data.shape), axis=-1)
    out_idx = np.nonzero(deln.find_simplex(idx) + 1)
    out_img = np.ones(seg_data.shape)
    out_img[out_idx] = 0

    # Save convex hull
    print("Saving convex hull...")
    out_hull = nib.Nifti1Image(out_img, seg_affine)
    nib.save(out_hull, out_convex_hull)


if __name__ == "__main__":
    create_convex_hull(
        binary_seg=snakemake.input.bin_seg,
        out_convex_hull=snakemake.output.convex_hull,
    )
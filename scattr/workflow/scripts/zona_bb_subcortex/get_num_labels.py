#!/usr/bin/env python
import nibabel as nib
import numpy as np


def get_num_labels(label_seg, out_num_labels):
    # Get number of labels
    img = nib.load(str(label_seg))
    img_data = img.get_fdata()
    num_labels = len(np.unique(img_data[img_data > 0]))

    # Write to file
    with open(out_num_labels, "w") as f:
        f.write(str(num_labels))


if __name__ == "__main__":
    get_num_labels(
        label_seg=snakemake.input.seg,  
        out_num_labels=snakemake.output.num_labels,  
    )

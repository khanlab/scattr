#!/usr/bin/env python
from nilearn.plotting import plot_roi

def show_segmented(t1w_path, seg_path, output_path_static, output_path_html, coords = None):
    plot_roi(seg_path, t1w_path, cut_coords = coords, output_file = output_path_static)
    view_img(seg_path, t1w_path, opacity=0.2).save_as_html(output_path_html)

if __name__ == "__main__":
    show_segmented(snakemake.input.t1w_image, snakemake.input.qc_labels, snakemake.output.segment_overlay, snakemake.output.segment_html)


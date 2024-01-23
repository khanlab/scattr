#!/usr/bin/env python
import matplotlib
from nilearn import plotting


def segmentation_qc(
    t1w_path, seg_path, output_path_static, output_path_html, wildcards
):
    matplotlib.use(backend="Agg")

    # Create HTML
    html_view = plotting.view_img(
        stat_map_img=seg_path,
        bg_img=t1w_path,
        opacity=0.3,
        cmap="viridis",
        dim=-1,
        symmetric_cmap=False,
        title="sub-{subject} dseg".format(**wildcards),
    )
    html_view.save_as_html(output_path_html)

    # Create static - mosaic with 7 cuts in each orientation
    display = plotting.plot_roi(
        roi_img=seg_path,
        bg_img=t1w_path,
        view_type="contours",
        display_mode="mosaic",
        cut_coords=7,
    )
    display.savefig(output_path_static)
    display.close()


if __name__ == "__main__":
    segmentation_qc(
        t1w_path=snakemake.input.t1w_image,
        seg_path=snakemake.input.qc_labels,
        output_path_static=snakemake.output.qc_png,
        output_path_html=snakemake.output.qc_html,
        wildcards=snakemake.wildcards,
    )

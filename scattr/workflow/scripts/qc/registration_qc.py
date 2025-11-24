#!/usr/bin/env python
import matplotlib
from nilearn.image import load_img
from nilearn.plotting import view_img
from niworkflows.viz.utils import compose_view, cuts_from_bbox, plot_registration


def registration_qc(
    moving_nii,
    fixed_nii,
    svg_cuts,
    out_svg,
    out_html,
    smk_wildcards,
):
    matplotlib.use(backend="Agg")

    # HTML output
    html_view = view_img(
        stat_map_img=moving_nii,
        bg_img=fixed_nii,
        opacity=0.5,
        cmap="viridis",
        dim=-1,
        symmetric_cmap="False",
        title="sub-{subject} registration".format(**smk_wildcards),
    )
    html_view.save_as_html(out_html)

    # SVG flicker output
    subj_img = load_img(fixed_nii)
    template_img = load_img(moving_nii)
    cuts = cuts_from_bbox(subj_img, svg_cuts)

    compose_view(
        plot_registration(  # Fixed image
            anat_nii=subj_img,
            div_id="fixed-img",
            estimate_brightness=True,
            cuts=cuts,
            label="Subject",
        ),
        plot_registration(
            anat_nii=template_img,
            div_id="mov-img",
            estimate_brightness=True,
            cuts=cuts,
            label="Template",
        ),
        out_file=out_svg,
    )


if __name__ == "__main__":
    registration_qc(
        moving_nii=snakemake.input.moving_nii,
        fixed_nii=snakemake.input.fixed_nii,
        svg_cuts=snakemake.params.cuts,
        out_svg=snakemake.output.qc_svg,
        out_html=snakemake.output.qc_html,
        smk_wildcards=snakemake.wildcards,
    )

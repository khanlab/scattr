#!/usr/bin/env python
import matplotlib
from nilearn import plotting


def registration_qc(
    template_t1w,
    native_t1w,
    out_png,
    out_html,
    smk_wildcards,
):
    matplotlib.use(backend="Agg")

    # HTML output
    html_view = plotting.view_img(
        stat_map_img=template_t1w,
        bg_img=native_t1w,
        opacity=0.5,
        cmap="viridis",
        dim=-1,
        symmetric_cmap="False",
        title="sub-{subject} registration".format(**smk_wildcards),
    )
    html_view.save_as_html(out_html)

    # PNG output
    display = plotting.plot_anat(native_t1w, display_mode="mosaic")
    display.add_edges(template_t1w, color="r")
    display.savefig(out_png)
    display.close()


if __name__ == "__main__":
    registration_qc(
        template_t1w=snakemake.input.template_t1w,  # noqa: F821
        native_t1w=snakemake.input.native_t1w,  # noqa: F821
        out_png=snakemake.output.qc_png,  # noqa: F821
        out_html=snakemake.output.qc_html,  # noqa: F821
        smk_wildcards=snakemake.wildcards,  # noqa: F821
    )

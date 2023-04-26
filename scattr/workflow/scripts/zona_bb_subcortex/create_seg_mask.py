#!/usr/env/bin python
import nibabel as nib


def create_seg_mask(dseg, out_mask):
    dseg_data = nib.load(dseg)
    affine, header = dseg_data.affine, dseg_data.header

    out_data = dseg_data.get_fdata() > 0
    out_nii = nib.Nifti1Image(dataobj=out_data, affine=affine, header=header)
    nib.save(out_nii, out_mask)


if __name__ == "__main__":
    create_seg_mask(
        dseg=snakemake.inputs.seg,  # noqa: F821
        out_mask=snakemake.output.mask,  # noqa: F821
    )

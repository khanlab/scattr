#!/usr/bin/env python
import nibabel as nib


def add_brainstem(dseg_mask, aparcaseg, out_mask):
    dseg_mask_data = nib.load(dseg_mask)
    aparcaseg_data = nib.load(aparcaseg)
    affine, header = aparcaseg_data.affine, aparcaseg_data.header

    out_data = (aparcaseg_data.get_fdata() == 16) | (dseg_mask_data.get_fdata() > 0)
    out_nii = nib.Nifti1Image(dataobj=out_data, affine=affine, header=header)
    nib.save(out_nii, out_mask)


if __name__ == "__main__":
    add_brainstem(
        dseg_mask=snakemake.input.mask,
        aparcaseg=snakemake.input.aparcaseg,
        out_mask=snakemake.output.mask,
    )

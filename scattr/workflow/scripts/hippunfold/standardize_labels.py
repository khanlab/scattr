		import glob, nibabel as nib
		from pathlib import Path
		
		sub = wildcards.subject
		base = Path(hippunfold_dir) / f"sub-{sub}" / "anat"
		patt_L = base / f"sub-{sub}_hemi-L_space-T1w_label-hipp_desc-subfields_atlas-multihist7_dseg.nii.gz"
		patt_R = base / f"sub-{sub}_hemi-R_space-T1w_label-hipp_desc-subfields_atlas-multihist7_dseg.nii.gz"
		
		hits_L = glob.glob(str(patt_L))
		hits_R = glob.glob(str(patt_R))
		
		if not hits_L or not hits_R:
			raise FileNotFoundError(f"Missing HippUnfold subfields:\n  {patt_L}\n  {patt_R}")
			
		imgL, imgR = nib.load(hits_L[0]), nib.load(hits_R[0])
		L = imgL.get_fdata().astype(int)
		R = imgR.get_fdata().astype(int)
		
		R[R > 0] += 100
		
		merged = L.copy()
		nz = R > 0
		merged[nz] = R[nz]
			
		Path(str(output.dseg)).parent.mkdir(parents=True, exist_ok=True)
		nib.Nifti1Image(merged, imgL.affine, imgL.header).to_filename(str(output.dseg))
		
		with open(str(log), "a") as LOG:
			LOG.write(f"Merged L={hits_L[0]}\n  R={hits_R[0]}\n  -> {output.dseg}\n")
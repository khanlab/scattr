hippunfold_dir = str(Path(config["output_dir"]) / "hippunfold")
if config.get("hippunfold_dir"):
	hippunfold_dir = config.get("hippunfold_dir")

log_dir = str(Path(config["output_dir"])/".logs"/"hippunfold")

bids_hippu_out = partial(
	bids,
	root=hippunfold_dir,
	datatype='anat',
	**inputs_t1w.subj_wildcards,
)


bids_hippu_log = partial(
	bids,
	root=log_dir,
	**inputs_t1w.subj_wildcards,
)

hippunfold_overlay_dir = (
	config.get("hippunfold_overlay_dir")
	or str(Path(config["output_dir"])/"hippunfold_bids")
)

bids_hippu_overlay = partial(
	bids,
	root=hippunfold_overlay_dir,
	datatype="anat",
	**inputs_t1w.subj_wildcards,
)

rule cp_hippu_tsv:
	input:
		hip_tsv=str(Path(workflow.basedir).parent / Path(config["hippunfold"]["tsv"]))
	output:
		hip_tsv=str(Path(hippunfold_overlay_dir) / "desc-HippUnfoldSubfields_dseg.tsv")
	threads: 1
	resources:
		mem_mb=2000,
		time=10,
	group:
		"dseg_tsv"
	shell:
		"""
		cp -v {input.hip_tsv} {output.hip_tsv} &> {log}
		"""

rule hippunfold_participant:
	input:
		bids_dir=lambda wc: config["bids_dir"]
	params:
		outdir=hippunfold_dir,
		modality=config["hippunfold"].get("modality", "T1w"),
		extra=config["hippunfold"].get("extra", "")
	output:
		sentinel=bids_hippu_log(suffix="hippunfold.done"),
	threads: config["hippunfold"].get("threads", 8)
	resources:
		mem_mb=16000,
		time=120,
	log:
		bids_hippu_log(suffix="hippunfold.log")
	group: "hippunfold"
	conda: config["hippunfold"].get("env", "envs/hippunfold.yml")
	shell:
		r"""
		hippunfold "{input.bids_dir}" "{params.outdir}" participant \
		--modality T1w \
		--cores all \
		--participant-label "{wildcards.subject}" \
		{params.extra} &> {log}
		
		touch {output.sentinel}
		"""


rule hippunfold_standardize_labels:
	input:
        	sentinel=rules.hippunfold_participant.output.sentinel,
	output:
		dseg=bids_hippu_out(
			space="T1w",
			desc="HippUnfoldSubfields",
			suffix="dseg.nii.gz",
		)
	threads: 4
	resources:
		mem_mb=16000, 
		time=20,
	log:
		bids_hippu_log(suffix="hippunfoldStandardizeLabels.log")
	conda:
		config["hippunfold"].get("env", "envs/hippunfold.yml")
	script:
		"../scripts/hippunfold/standardize_labels.py"


rule hippunfold_export_overlay:
	input:
		dseg=rules.hippunfold_standardize_labels.output.dseg	
	output:
		overlay=bids_hippu_overlay(
			space="T1w",
			desc="HippUnfoldSubfields",
			suffix="dseg.nii.gz",
		)
	threads: 1
	resources:
		mem_mb=2000,
		time=5,
	log:
		bids_hippu_log(suffix="exportOverlay.log")
	shell:
		r"""
		mkdir -p "$(dirname {output.overlay})"
		cp -v {input.dseg} {output.overlay} &> {log}
		"""

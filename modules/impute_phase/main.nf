//Processes
process impute_phase {
	tag "${run_id}_${sample_id}_${chr_id}"
	publishDir "${params.outdir}/${run_id}_${sample_id}/imputation/phase", mode:"copy"

	input:
	tuple val(chr_id), file(chunk_txt), file(ref_bin), val(run_id), val(sample_id), file(sample_cram), file(sample_index), file(ref_fasta), file(fasta_index)

	output:
	tuple val(chr_id), val(run_id), val(sample_id), path("*.bcf"), path("*.bcf.csi"), emit: phase_out

	script:
	"""
	python /wrapper/GLIMPSE_phase.py -S ${run_id}_${sample_id} -C ${chr_id} -I ${sample_cram} -R ${ref_bin} -c ${chunk_txt} -F ${ref_fasta}
	"""
}

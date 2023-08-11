process impute_concat {
	tag "${run_id}_${sample_id}"
	publishDir "${params.outdir}/${run_id}_${sample_id}/imputation/final", mode:"copy"

	input:
	tuple val(run_id), val(sample_id), path(bcf), path(index)

	output:
	tuple val(run_id), val(sample_id), path("*.final.bcf"), path("*.final.bcf.csi"), emit: concat_out

	script:
	"""
	bcftools concat -a ${bcf} -o ${run_id}_${sample_id}.final.bcf
	bcftools index -f ${run_id}_${sample_id}.final.bcf
	"""
}

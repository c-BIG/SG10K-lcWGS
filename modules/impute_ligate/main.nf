process impute_ligate {
	tag "${run_id}_${sample_id}_${chr_id}"
	publishDir "${params.outdir}/${run_id}_${sample_id}/imputation/ligate", mode:"copy"

	input:
	tuple val(chr_id), val(run_id), val(sample_id), path(bcf), path(index)

	output:
	tuple val(run_id), val(sample_id), path("*.merged.bcf"), path("*.merged.bcf.csi"), emit:ligate_out
	path ("*.txt"), emit:txt_out

	script:
	"""
	for f in ${bcf}; do echo "\${f}"; done | sort -V > ${run_id}_${sample_id}_${chr_id}.txt
	GLIMPSE2_ligate --input ${run_id}_${sample_id}_${chr_id}.txt --output ${run_id}_${sample_id}_${chr_id}.merged.bcf
	bcftools index -f ${run_id}_${sample_id}_${chr_id}.merged.bcf
	"""
}

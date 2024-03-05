process ligatecommon {
          debug true

	tag ${chr_no}
	publishDir "${params.outdir}/${chr_no}/ligate", mode:"copy"

	input:
	tuple val(chr_no), path(bcf), path(index)

	output:
	// tuple val(chr_no), path("*.merged.bcf"), path("*.merged.bcf.csi"), emit:ligate_out
	path ("${chr_no}/${chr_no}_phase_common.txt"), emit:txt_out

	script:
	"""
        ls -1v ${chr_no}/SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_*.common.bcf > ${chr_no}/${chr_no}_phase_common.txt

  """
}

}

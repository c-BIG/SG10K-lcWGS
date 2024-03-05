process ligatecommon {

	tag {chr_no}

	input:
	tuple val(chr_no), path(bcf), path(index)

	output:
	path("${chr_no}/${chr_no}_phase_common.txt"), emit:common_list

	script:
	"""
        ls -1v ${chr_no}/SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_*.common.bcf > ${chr_no}/${chr_no}_phase_common.txt

        """
}

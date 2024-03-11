process testconcat {

	tag {chr_no}

	input:
	//tuple val(chr_no), path(bcf), path(index)
        tuple val(chr_no), path(bcf_files)

	output:
	tuple val(chr_no), path("${chr_no}/10k.${chr_no}.rare.full.bcf*")
        //tuple val(chr_no), path("${chr_no}/${chr_no}_phase_rare.txt"), emit:common_list

	script:
	"""
        mkdir -p ${chr_no}

        OUT=10k.${chr_no}.rare.full.bcf

        #ls -1v ${chr_no}/SG10K_Health_r5.3.2.n6686.shapeit5.fil`tered.${chr_no}.chunk_*.common.bcf > ${chr_no}/${chr_no}_phase_common.txt
        ls -1v *${chr_no}.chunk_*.rare.bcf > ${chr_no}/${chr_no}_phase_rare.txt
        echo TEST BCF ligate ${chr_no} ${bcf_files} ${chr_no}/${chr_no}_phase_rare.txt > ${chr_no}/\${OUT}
        echo TEST BCF ligate ${chr_no} ${bcf_files} > ${chr_no}/\${OUT}.csi

        """
}

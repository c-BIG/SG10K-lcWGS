process phasing {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf), path(tbi), path(chunks), path(gmap)

    output:
    tuple val(chr_no), path("${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_common.bcf"), path("${chr_no}/${chr_no}_phase_common.txt"), path("${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.shapeit5_common_ligate.bcf")
    // tuple val(chr_no), path("SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_common.bcf"), emit: phased_bcf
    // tuple val(chr_no), path("${chr_no}_phase_common.txt"), emit: phased_bcf_list
    // tuple val(chr_no), path(SG10K_Health_r5.3.2.n9770.${chr_no}.shapeit5_common_ligate.bcf), emit: ligate_bcf

    script:

    bcf_threads = (int) Math.ceil(task.cpus*2)
    """

    #phasing per chr by chunks

    while read LINE;
    do
    REG=\$(echo \$LINE | awk '{ print \$3; }')
    CHUNK_NBR=\$(echo \$LINE | awk '{ print \$1; }')
    OUT=SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_\${CHUNK_NBR}.shapeit5_common.bcf

    phase_common \
        --input ${vcf} \
        --map ${gmap} \
        --output ${chr_no}/\${OUT} \
        --thread ${bcf_threads} \
        --filter-maf 0.001 \
        --region \${REG} && bcftools index -f ${chr_no}/\${OUT} --threads ${bcf_threads} \
    done < ${chunks}

ls -1v ${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_common.bcf > ${chr_no}/${chr_no}_phase_common.txt

    ligate \
    --input ${chr_no}/${chr_no}_phase_common.txt \
    --output ${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.shapeit5_common_ligate.bcf \
    --thread ${bcf_threads} \
    --index
    """
}

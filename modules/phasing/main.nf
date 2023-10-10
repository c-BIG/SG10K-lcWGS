process phasing {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf), path(tbi), path(chunks), path(gmap)

    output:
    tuple val(chr_no), path("${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_common.bcf"), path("${chr_no}/${chr_no}_phase_common.txt"), path("${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.shapeit5_common_ligate.bcf"), path("${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_rare.bcf"), path("${chr_no}/SG10K_Health_r5.3.2.n9770_${chr_no}.full.shapeit5_rare.bcf")

    // tuple val(chr_no), path("SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_common.bcf"), emit: phased_bcf
    // tuple val(chr_no), path("${chr_no}_phase_common.txt"), emit: phased_bcf_list
    // tuple val(chr_no), path(SG10K_Health_r5.3.2.n9770.${chr_no}.shapeit5_common_ligate.bcf), emit: ligate_bcf

    script:

    bcf_threads = (int) Math.ceil(task.cpus*2)
    """

    #phasing common variants per chr by chunks

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
        --region \${REG} && bcftools index -f ${chr_no}/\${OUT} --threads ${bcf_threads}
    done < ${chunks}

    # list phased common varinats files per chr

    ls -1v ${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_common.bcf > ${chr_no}/${chr_no}_phase_common.txt

    # ligate phased common variants per chr

    ligate \
    --input ${chr_no}/${chr_no}_phase_common.txt \
    --output ${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.shapeit5_common_ligate.bcf \
    --thread ${bcf_threads} \
    --index

    # phasing rare variants per chr by chunks

    SCAF=${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.shapeit5_common_ligate.bcf

    while read LINE;
    do
    CHUNK_NBR=\$(echo \$LINE | awk '{ print \$1; }')

    SCAFFOLD_REG=\$(echo \$LINE | awk '{ print \$3; }')
    SCAFFOLD_REG_START=\$(echo \${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f1)
    SCAFFOLD_REG_END=\$(echo \${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f2)
    SCAFFOLD_REG_NAME=${chr_no}_\${SCAFFOLD_REG_START}_\${SCAFFOLD_REG_END}
			
    INPUT_REG=\$(echo \$LINE | awk '{ print \$4; }')
    INPUT_REG_START=\$(echo \${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
    INPUT_REG_END=\$(echo \${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
    INPUT_REG_NAME=${chr_no}_\${INPUT_REG_START}_\${INPUT_REG_END}

    OUT=SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_\${CHUNK_NBR}.shapeit5_rare.bcf

    phase_rare \
    --input-plain ${vcf} \
    --map ${gmap} \
    --output ${chr_no}/\$OUT \
    --thread ${bcf_threads} \
    --scaffold \$SCAF \
    --scaffold-region \$SCAFFOLD_REG \
    --input-region \$INPUT_REG && bcftools index -f \$OUT --threads ${bcf_threads}
    done < ${chunks}

    # list chr wide phased data

    ls -1v ${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_rare.bcf > ${chr_no}/concat_list_${chr_no}.txt

    # concatenate chunks

    OUT_CONCAT=${chr_no}/SG10K_Health_r5.3.2.n9770_${chr_no}.full.shapeit5_rare.bcf

    bcftools concat -n -f ${chr_no}/concat_list_${chr_no}.txt -o ${OUT_CONCAT} && bcftools index ${OUT_CONCAT} --threads ${bcf_threads}

    """
}

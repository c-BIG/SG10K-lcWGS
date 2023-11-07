process phasing {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf), path(tbi), path(chunks), path(gmap), path(scaffold), path(scaffold_idx)

    output:
    tuple val(chr_no), path("${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_rare.bcf"), path("${chr_no}/concat_list_${chr_no}.txt")


    script:

    """

    mkdir -p ${chr_no}

    # phasing rare variants per chr by chunks

    ## SCAF=${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.shapeit5_common_ligate.bcf

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
        --input ${vcf} \
        --map ${gmap} \
        --output ${chr_no}/\${OUT} \
        --thread ${task.cpus} \
        --scaffold ${scaffold} \
        --scaffold-region \${SCAFFOLD_REG} \
        --input-region \${INPUT_REG} && bcftools index -f ${chr_no}/\${OUT} --threads ${task.cpus}
    done < ${chunks}

    # list chr wide phased data

    ls -1v ${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_rare.bcf > ${chr_no}/concat_list_${chr_no}.txt

    """
}

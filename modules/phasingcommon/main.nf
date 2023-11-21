process phasingcommon {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf), path(tbi), path(chunks), path(gmap)

    output:
    tuple val(chr_no), path("${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_common.bcf*"), path("${chr_no}/${chr_no}_phase_common.txt") 


    script:

    """

    mkdir -p ${chr_no}


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
        --thread ${task.cpus} \
        --filter-maf 0.001 \
        --region \${REG} && bcftools index -f ${chr_no}/\${OUT} --threads ${task.cpus}
    done < ${chunkscommon}

    # list phased common varinats files per chr

    ls -1v ${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.chunk_*.shapeit5_common.bcf > ${chr_no}/${chr_no}_phase_common.txt

    """
}

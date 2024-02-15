process phasingcommon {

    tag "${chunk_no}_${chr_no}_${chunk_region}"

    input:
    tuple val(chunk_no), val(chr_no), val(chunk_region), path(vcf), path(tbi), path(gmap)

    output:
    tuple val(chunk_no), val(chr_no), val(chunk_region), path("${chr_no}/SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_*.common.bcf*")


    script:

    """
    mkdir -p ${chr_no}

    #phasing common variants per chr by chunks

    OUT=SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_${chunk_no}.common.bcf

    phase_common \
        --input ${vcf} \
        --map ${gmap} \
        --output ${chr_no}/\${OUT} \
        --thread ${task.cpus} \
        --filter-maf 0.001 \
        --region ${chunk_region} && bcftools index -f ${chr_no}/\${OUT} --threads ${task.cpus}

    """
}

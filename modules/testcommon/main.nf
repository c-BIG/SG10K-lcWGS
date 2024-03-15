process testcommon {

    debug true
    tag "${chr_no}_${chunk_no}_${chunk_region}"

    input:
    tuple val(chr_no), val(chunk_no), val(chunk_region), path(bcf), path(csi), path(gmap)
    val(PREFIX)

    output:
    tuple val(chr_no), path("${chr_no}/${PREFIX}.${chr_no}.chunk_*.common.bcf*"), emit: common_bcf_files

    script:

    """
    mkdir -p ${chr_no}

    #phasing common variants per chr and chunks


    phase_common \
        --input ${bcf} \
        --map ${gmap} \
        --output ${chr_no}/${PREFIX}.${chr_no}.chunk_${chunk_no}.common.bcf \
        --thread ${task.cpus} \
        --filter-maf 0.001 \
        --region ${chunk_region} && bcftools index -f ${chr_no}/${PREFIX}.${chr_no}.chunk_${chunk_no}.common.bcf --threads ${task.cpus}
    """
}

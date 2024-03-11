process testcommon {

    debug true
    tag "${chr_no}_${chunk_no}_${chunk_region}"

    input:
    tuple val(chr_no), val(chunk_no), val(chunk_region), path(bcf), path(csi), path(gmap)

    output:
    // tuple val(chr_no), val(chunk_no), val(chunk_region), path("${chr_no}/SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_*.common.bcf*")
    // tuple val(chr_no), path("${chr_no}/SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_*.common.bcf*")
    tuple val(chr_no), path("${chr_no}/10k.${chr_no}.chunk_*.common.bcf*")

    script:

    """
    mkdir -p ${chr_no}

    #phasing common variants per chr by chunks

    #OUT=SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_${chunk_no}.common.bcf
    #OUT_idx=SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_${chunk_no}.common.bcf.csi
    OUT=10k.${chr_no}.chunk_${chunk_no}.common.bcf
    OUT_idx=10k.${chr_no}.chunk_${chunk_no}.common.bcf.csi

    echo TEST BCF ${chunk_no} ${chunk_region} ${bcf} > ${chr_no}/\${OUT}
    echo TEST IDX ${chunk_no} ${chunk_region} ${bcf} > ${chr_no}/\${OUT_idx}

    """
}

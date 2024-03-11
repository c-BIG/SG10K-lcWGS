process testrare {

    debug true
    tag "${chr_no}_${chunk_no}_${scaffold_reg}_${input_reg}"

    input:
    tuple val(chr_no), val(chunk_no), val(scaffold_reg), val(input_reg), path(bcf), path(csi), path(gmap), path(scaffold)

    output:
    // tuple val(chr_no), val(chunk_no), val(chunk_region), path("${chr_no}/SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_*.common.bcf*")
    // tuple val(chr_no), path("${chr_no}/SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_*.common.bcf*")
    tuple val(chr_no), path("${chr_no}/10k.${chr_no}.chunk_*.rare.bcf*")
    //tuple val(chr_no), path("${chr_no}/test-rare.txt")

    script:

    """
    mkdir -p ${chr_no}

    #phasing common variants per chr by chunks

    #OUT=SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_${chunk_no}.common.bcf
    #OUT_idx=SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_${chunk_no}.common.bcf.csi
    OUT=10k.${chr_no}.chunk_${chunk_no}.rare.bcf
    OUT_idx=10k.${chr_no}.chunk_${chunk_no}.rare.bcf.csi

    echo ${chr_no}, ${chunk_no}, ${scaffold_reg}, ${input_reg}, ${bcf}, ${csi}, ${gmap}, ${scaffold[0]} >${chr_no}/\${OUT}
    # echo ${chr_no}, ${chunk_no}, ${scaffold_reg}, ${input_reg}, ${bcf}, ${csi}, ${gmap}, ${scaffold[0]} >${chr_no}/\${OUT_idx}
    #echo TEST BCF --input ${bcf} ${csi} --map ${gmap} --scaffold ${scaffold} --scaffold-region ${scaffold_reg} --input-region ${input_reg} ${chunk_no}  > ${chr_no}/\${OUT}
    #echo TEST IDX ${chr_no} ${chunk_no} ${scaffold_reg} ${input_reg} ${bcf} > ${chr_no}/\${OUT_idx}

    """
}

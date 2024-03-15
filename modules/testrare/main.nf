process testrare {

    debug true
    tag "${chr_no}_${chunk_no}_${scaffold_reg}_${input_reg}"

    input:
    tuple val(chr_no), val(chunk_no), val(scaffold_reg), val(input_reg), path(bcf), path(csi), path(gmap), path(scaffold)
    val(PREFIX)

    output:
    tuple val(chr_no), path("${chr_no}/${PREFIX}.${chr_no}.chunk_*.rare.bcf*"), emit: rare_scaffold_files
    //tuple val(chr_no), path("${chr_no}/test-rare.txt")

    script:

    """
    mkdir -p ${chr_no}

    #phasing rare variants per chr by chunks

    OUT=${PREFIX}.${chr_no}.chunk_${chunk_no}.rare.bcf

    phase_rare \
        --input ${bcf} \
        --map ${gmap} \
        --output ${chr_no}/\${OUT} \
        --thread ${task.cpus} \
        --scaffold ${scaffold[0]} \
        --scaffold-region ${scaffold_reg} \
        --input-region ${input_reg} && bcftools index -f ${chr_no}/\${OUT} --threads ${task.cpus}
    """
}

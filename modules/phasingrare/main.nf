process phasingrare {

    tag "${chr_no}_${chunk_no}_${scaffold_region}_${input_region}"

    input:
    tuple val(chr_no), val(chunk_no), val(scaffold_region), val(input_region)
    path(vcf)
    path(tbi)
    path(gmap)
    path(scaffold)
    path(scaffold_idx)

    output:
    tuple val(chr_no), path("${chr_no}/SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_*.common.bcf*")


    script:

    """

    mkdir -p ${chr_no}

    # phasing rare variants per chr by chunks

    ## SCAF=${chr_no}/SG10K_Health_r5.3.2.n9770.${chr_no}.shapeit5_common_ligate.bcf

    OUT=SG10K_Health_r5.3.2.n6686.shapeit5.filtered.${chr_no}.chunk_${chunk_no}.rare.bcf

    phase_rare \
        --input ${vcf} \
        --map ${gmap} \
        --output ${chr_no}/\${OUT} \
        --thread ${task.cpus} \
        --scaffold ${scaffold} \
        --scaffold-region ${scaffold_region} \
        --input-region ${input_region} && bcftools index -f ${chr_no}/\${OUT} --threads ${task.cpus}

    """
}

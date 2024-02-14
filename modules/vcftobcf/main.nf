process vcftobcf {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf), path(tbi)

    output:
    tuple val(chr_no), path("*.bcf"), path("*.bcf.csi"), emit: bcf_out


    script:

    """

    #convert to bcf of qced SG10K Health 5.3.2 vcf which exported from hail metrics table

    bcftools view -Ob --threads ${task.cpus*2} -o SG10K_Health_r5.3.shapeit5-qc-filtered.${chr_no}.bcf ${vcf}
    bcftools index -f SG10K_Health_r5.3.shapeit5-qc-filtered.${chr_no}.bcf --threads ${task.cpus*2}

    """
}

process refpanel {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf)
    path reheader

    output:
    tuple val(chr_no), path("sg10k_chr22.reheader.vcf.gz")

    script:

    bcf_threads = (int) Math.ceil(task.cpus*2)
    """

    bcftools reheader -h ${reheader} ${vcf} --threads ${bcf_threads} -o sg10k_chr22.reheader.vcf
    bgzip sg10k_chr22.reheader.vcf

    """
}

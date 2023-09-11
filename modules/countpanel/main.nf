process countpanel {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf_sg10k), path(tbi_sg10k), path(vcf_1kg), path(tbi_1kg)

    output:
    tuple val(chr_no), path("${chr_no}.sg10k.count.txt"), path("${chr_no}.1kg.count.txt"), path("${chr_no}_isec/*")

    script:

    """

    #create number of variants count in sg10k and 1kg reference panel
    #create the number of common variants in sg10k and 1kg reference panel
    bcftools view -r ${chr_no} ${vcf_sg10k} --threads ${task.cpus} | grep -v -c '^#' >${chr_no}.sg10k.count.txt &
    bcftools view -r ${chr_no} ${vcf_1kg} --threads ${task.cpus} | grep -v -c '^#' >${chr_no}.1kg.count.txt &
    bcftools isec -p ${chr_no}_isec -n=2 -w1 ${vcf_sg10k} ${vcf_1kg}
    wc -l ${chr_no}_isec/sites.txt >${chr_no}_isec/${chr_no}.isec.count.txt
    """
}

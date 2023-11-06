process af_dist {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf_sg10k), path(tbi_sg10k), path(vcf_1kg), path(tbi_1kg)

    output:
    tuple val(chr_no), path("${1kg_AF_MAF_bin_out}"), path("${sg10k_AF_MAF_bin_out}")

    script:

    """

    touch 1kg_${chr_no}.AF.MAF_bin_count_all.txt
    touch sg10k_${chr_no}.MAF_bin_count_all.txt
    touch 1kg_${chr_no}.AF.MAF_bin_count_output.txt
    touch sg10k_${chr_no}.MAF_bin_count_ouput.txt

    1kg_AF_MAF_bin="1kg_${chr_no}.AF.MAF_bin_count_all.txt"
    sg10k_AF_MAF_bin="sg10k_${chr_no}.MAF_bin_count_all.txt"
    1kg_AF_MAF_bin_out="1kg_${chr_no}.AF.MAF_bin_count_output.txt"
    sg10k_AF_MAF_bin_out="sg10k_${chr_no}.MAF_bin_count_ouput.txt"


    #create number of variants count in sg10k and 1kg reference panel
    #bin the MAF
    bcftools view -H ${vcf_sg10k} --threads ${task.cpus} |cut -f8 |cut -f2 -d';'|sed 's/AF=//g' >${sg10k_AF_MAF_bin}

    bcftools +fill-tags ${vcf_1kg} -- -t AF | bcftools view --threads ${task.cpus} --drop-genotypes -Oz -o 1kg_${chr_no}_AF.vcf.gz -
    tabix 1kg_${chr_no}_AF.vcf.gz

    bcftools view -H ${vcf_sg10k} --threads ${task.cpus} |cut -f8 |cut -f2 -d';'|sed 's/AF=//g' >${sg10k_AF_MAF_bin}
    bcftools view -H 1kg_${chr_no}_AF.vcf.gz --threads ${task.cpus} |cut -f8|sed 's/AF=//g' >${1kg_AF_MAF_bin}

    awk '{ if ( \$1 <= 0.00000 ) print \$0 }' ${1kg_AF_MAF_bin}  >${chr_no}_1kg.AF.MAF_bin_count_0-00000.txt
    awk '{ if ( \$1 > 0.00000 && \$1 <= 0.00100) print \$0 }' ${1kg_AF_MAF_bin} >${chr_no}_1kg.AF.MAF_bin_count_0-00100.txt
    awk '{ if ( \$1 > 0.00100 && \$1 <= 0.00200) print \$0 }' ${1kg_AF_MAF_bin} >${chr_no}_1kg.AF.MAF_bin_count_0-00200.txt
    awk '{ if ( \$1 > 0.00200 && \$1 <= 0.00500) print \$0 }' ${1kg_AF_MAF_bin} >${chr_no}_1kg.AF.MAF_bin_count_0-00500.txt
    awk '{ if ( \$1 > 0.00500 && \$1 <= 0.01000) print \$0 }' ${1kg_AF_MAF_bin} >${chr_no}_1kg.AF.MAF_bin_count_0-01000.txt
    awk '{ if ( \$1 > 0.01000 && \$1 <= 0.05000) print \$0 }' ${1kg_AF_MAF_bin} >${chr_no}_1kg.AF.MAF_bin_count_0-05000.txt
    awk '{ if ( \$1 > 0.05000 && \$1 <= 0.10000) print \$0 }' ${1kg_AF_MAF_bin} >${chr_no}_1kg.AF.MAF_bin_count_0-10000.txt
    awk '{ if ( \$1 > 0.10000 && \$1 <= 0.20000) print \$0 }' ${1kg_AF_MAF_bin} >${chr_no}_1kg.AF.MAF_bin_count_0-20000.txt
    awk '{ if ( \$1 > 0.20000 && \$1 <= 0.50000) print \$0 }' ${1kg_AF_MAF_bin} >${chr_no}_1kg.AF.MAF_bin_count_0-50000.txt


    awk '{ if ( \$1 <= 0.00000 ) print \$0 }' ${sg10k_AF_MAF_bin}  >${chr_no}_sg10k.AF.MAF_bin_count_0-00000.txt

    cat ${chr_no}_1kg.AF.MAF_bin_count*.txt >${1kg_AF_MAF_bin_out}
    cat ${chr_no}_sg10k.AF.MAF_bin_count*.txt >${sg10k_AF_MAF_bin_out}

    """
}

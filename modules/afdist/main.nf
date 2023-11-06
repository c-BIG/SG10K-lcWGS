process af_dist {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf_sg10k), path(tbi_sg10k), path(vcf_1kg), path(tbi_1kg)

    output:
    # tuple val(chr_no), path("1kg_${chr_no}.AF.MAF_bin_count_output.txt"), path("sg10k_${chr_no}.MAF_bin_count_ouput.txt")
    tuple val(chr_no), path("1kg_${chr_no}.AF.MAF_bin_count_output.txt")

    script:

    """

    #create number of variants count in sg10k and 1kg reference panel
    #bin the MAF
    bcftools +fill-tags ${vcf_1kg} -- -t AF | bcftools view --threads ${task.cpus} --drop-genotypes -Oz -o 1kg_${chr_no}_AF.vcf.gz -
    tabix 1kg_${chr_no}_AF.vcf.gz

    #bcftools view -H ${vcf_sg10k} --threads ${task.cpus} |cut -f8 |cut -f2 -d';'|sed 's/AF=//g' >sg10k_${chr_no}.MAF_bin_count_all.txt
    bcftools view -H 1kg_${chr_no}_AF.vcf.gz --threads ${task.cpus} |cut -f8|sed 's/AF=//g' >1kg_${chr_no}.AF.MAF_bin_count_all.txt

    awk '{ if ( \$1 <= 0.00000 ) print \$0 }' 1kg_${chr_no}.AF.MAF_bin_count_all.txt  >${chr_no}_1kg.AF.MAF_bin_count_0-00000.txt
    awk '{ if ( \$1 > 0.00000 && \$1 <= 0.00100) print \$0 }' 1kg_${chr_no}.AF.MAF_bin_count_all.txt >${chr_no}_1kg.AF.MAF_bin_count_0-00100.txt
    awk '{ if ( \$1 > 0.00100 && \$1 <= 0.00200) print \$0 }' 1kg_${chr_no}.AF.MAF_bin_count_all.txt >${chr_no}_1kg.AF.MAF_bin_count_0-00200.txt
    awk '{ if ( \$1 > 0.00200 && \$1 <= 0.00500) print \$0 }' 1kg_${chr_no}.AF.MAF_bin_count_all.txt >${chr_no}_1kg.AF.MAF_bin_count_0-00500.txt
    awk '{ if ( \$1 > 0.00500 && \$1 <= 0.01000) print \$0 }' 1kg_${chr_no}.AF.MAF_bin_count_all.txt >${chr_no}_1kg.AF.MAF_bin_count_0-01000.txt
    awk '{ if ( \$1 > 0.01000 && \$1 <= 0.05000) print \$0 }' 1kg_${chr_no}.AF.MAF_bin_count_all.txt >${chr_no}_1kg.AF.MAF_bin_count_0-05000.txt
    awk '{ if ( \$1 > 0.05000 && \$1 <= 0.10000) print \$0 }' 1kg_${chr_no}.AF.MAF_bin_count_all.txt >${chr_no}_1kg.AF.MAF_bin_count_0-10000.txt
    awk '{ if ( \$1 > 0.10000 && \$1 <= 0.20000) print \$0 }' 1kg_${chr_no}.AF.MAF_bin_count_all.txt >${chr_no}_1kg.AF.MAF_bin_count_0-20000.txt
    awk '{ if ( \$1 > 0.20000 && \$1 <= 0.50000) print \$0 }' 1kg_${chr_no}.AF.MAF_bin_count_all.txt >${chr_no}_1kg.AF.MAF_bin_count_0-50000.txt


    #awk '{ if ( \$1 <= 0.00000 ) print \$0 }' sg10k_${chr_no}.MAF_bin_count_all.txt  >${chr_no}_sg10k.AF.MAF_bin_count_0-00000.txt

    cat ${chr_no}_1kg.AF.MAF_bin_count*.txt >1kg_${chr_no}.AF.MAF_bin_count_output.txt
    #cat ${chr_no}_sg10k.AF.MAF_bin_count*.txt >sg10k_${chr_no}.MAF_bin_count_ouput.txt

    """
}

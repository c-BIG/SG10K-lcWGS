process af_dist {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf_sg10k), path(tbi_sg10k), path(vcf_1kg), path(tbi_1kg)

    output:
    tuple val(chr_no), path("1kg_${chr_no}.AF.MAF_bin_count_output.txt"), path("1kg_${chr_no}.AF.MAF_count_all.txt"), path("sg10k_${chr_no}.MAF_count_all.txt"), path("sg10k_${chr_no}.MAF_bin_count_ouput.txt")

    script:

    """

    #create number of variants count in sg10k and 1kg reference panel
    #bin the MAF
    bcftools +fill-tags ${vcf_1kg} -- -t AF | bcftools view --threads ${task.cpus} --drop-genotypes -Oz -o 1kg_${chr_no}_AF.vcf.gz -
    tabix 1kg_${chr_no}_AF.vcf.gz

    bcftools view -H ${vcf_sg10k} --threads ${task.cpus} |cut -f8 |cut -f2 -d';'|sed 's/AF=//g' >sg10k_${chr_no}.AF.MAF_count_all.txt
    bcftools view -H 1kg_${chr_no}_AF.vcf.gz --threads ${task.cpus} |cut -f8|sed 's/AF=//g' >1kg_${chr_no}.AF.MAF_count_all.txt

    awk '{ if ( \$1 <= 0.00000 ) print \$0 }'                 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_00000.txt
    awk '{ if ( \$1 > 0.00000 && \$1 <= 0.00020) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_0002.txt
    awk '{ if ( \$1 > 0.00020 && \$1 <= 0.00050) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_0005.txt
    awk '{ if ( \$1 > 0.00050 && \$1 <= 0.00100) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_001.txt
    awk '{ if ( \$1 > 0.00100 && \$1 <= 0.00200) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_002.txt
    awk '{ if ( \$1 > 0.00200 && \$1 <= 0.00500) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_005.txt
    awk '{ if ( \$1 > 0.00500 && \$1 <= 0.01000) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_01.txt
    awk '{ if ( \$1 > 0.01000 && \$1 <= 0.02000) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_02.txt
    awk '{ if ( \$1 > 0.02000 && \$1 <= 0.05000) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_05.txt
    awk '{ if ( \$1 > 0.05000 && \$1 <= 0.10000) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_1.txt
    awk '{ if ( \$1 > 0.10000 && \$1 <= 0.15000) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_15.txt
    awk '{ if ( \$1 > 0.15000 && \$1 <= 0.20000) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_2.txt
    awk '{ if ( \$1 > 0.20000 && \$1 <= 0.30000) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_3.txt
    awk '{ if ( \$1 > 0.30000 && \$1 <= 0.40000) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_4.txt
    awk '{ if ( \$1 > 0.40000 && \$1 <= 0.50000) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_5.txt
    awk '{ if ( \$1 > 0.50000 && \$1 <= 0.60000) print \$0 }' 1kg_{chr_no}.AF.MAF_count_all.txt  |wc -l >1kg_{chr_no}.AF.MAF_bin_count_6.txt

    awk '{ if ( \$1 <= 0.00000 ) print \$0 }'                 sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_00000.txt
    awk '{ if ( \$1 > 0.00000 && \$1 <= 0.00020) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_0002.txt
    awk '{ if ( \$1 > 0.00020 && \$1 <= 0.00050) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_0005.txt
    awk '{ if ( \$1 > 0.00050 && \$1 <= 0.00100) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_001.txt
    awk '{ if ( \$1 > 0.00100 && \$1 <= 0.00200) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_002.txt
    awk '{ if ( \$1 > 0.00200 && \$1 <= 0.00500) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_005.txt
    awk '{ if ( \$1 > 0.00500 && \$1 <= 0.01000) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_01.txt
    awk '{ if ( \$1 > 0.01000 && \$1 <= 0.02000) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_02.txt
    awk '{ if ( \$1 > 0.02000 && \$1 <= 0.05000) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_05.txt
    awk '{ if ( \$1 > 0.05000 && \$1 <= 0.10000) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_1.txt
    awk '{ if ( \$1 > 0.10000 && \$1 <= 0.15000) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_15.txt
    awk '{ if ( \$1 > 0.15000 && \$1 <= 0.20000) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_2.txt
    awk '{ if ( \$1 > 0.20000 && \$1 <= 0.30000) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_3.txt
    awk '{ if ( \$1 > 0.30000 && \$1 <= 0.40000) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_4.txt
    awk '{ if ( \$1 > 0.40000 && \$1 <= 0.50000) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_5.txt
    awk '{ if ( \$1 > 0.50000 && \$1 <= 0.60000) print \$0 }' sg10k_{chr_no}.AF.MAF_count_all.txt  |wc -l >sg10k_{chr_no}.AF.MAF_bin_count_6.txt

	
    cat ${chr_no}_1kg.AF.MAF_bin_count*.txt >1kg_${chr_no}.AF.MAF_bin_count_output.txt
    cat ${chr_no}_sg10k.AF.MAF_bin_count*.txt >sg10k_${chr_no}.MAF_bin_count_ouput.txt

    """
}

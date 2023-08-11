process refpanel {

    tag { chr_no }
    publishDir "${params.publishdir}/ref_panel", mode:"copy"

    input:
    tuple val(chr_no), path(vcf), path(tbi), path(gmap)
    val(BIN)

    output:
    tuple val(chr_no), path("binary/${chr_no}_bin/${BIN}*"), path("chunks/${chr_no}_all.txt")

    script:

    """

    mkdir -p normalized
    mkdir -p tagged
    mkdir -p nogt
    mkdir -p chunks
    mkdir -p binary/${chr_no}_bin

    bcftools norm -m -any ${vcf} -Ou --threads 8 | \
    bcftools view -m 2 -M 2 -v snps,indels --threads 4 -Ob -o normalized/${BIN}_${chr_no}.normalized.bcf
    bcftools index -f normalized/${BIN}_${chr_no}.normalized.bcf --threads 8

    #adding allele number and allele count
    bcftools +fill-tags normalized/${BIN}_${chr_no}.normalized.bcf -Ob -o tagged/${BIN}_${chr_no}.tagged.bcf -- -t AN,AC
    bcftools index -f tagged/${BIN}_${chr_no}.tagged.bcf --threads 8

    #create bcf with no GT infor for chunking
    bcftools view -G -Ob -o nogt/${BIN}_${chr_no}.nogt.bcf tagged/${BIN}_${chr_no}.tagged.bcf
    bcftools index -f nogt/${BIN}_${chr_no}.nogt.bcf --threads 8

    #Chunking reference panel
    GLIMPSE2_chunk \
        --input nogt/${BIN}_${chr_no}.nogt.bcf \
        --region ${chr_no} \
        --output chunks/"${chr_no}_all.txt" \
        --map ${gmap} \
        --sequential

    #split reference panel into bin
    while IFS="" read -r LINE || [ -n "\$LINE" ];
    do
    printf -v ID "%02d" \$(echo \$LINE | cut -d" " -f1)
    IRG=\$(echo \$LINE | cut -d" " -f3)
    ORG=\$(echo \$LINE | cut -d" " -f4)

    GLIMPSE2_split_reference \
        --reference tagged/${BIN}_${chr_no}.tagged.bcf \
        --map ${gmap} \
        --input-region \${IRG} \
        --output-region \${ORG} \
        --output binary/${chr_no}_bin/${BIN}
    done < chunks/"${chr_no}_all.txt"

    """
}

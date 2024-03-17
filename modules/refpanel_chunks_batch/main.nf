process refpanel_chunks_batch {

    tag { chr_no }

    input:
    tuple val(chr_no), path(vcf), path(tbi), path(gmap)
    val BIN
    path reheader

    output:
    tuple val(chr_no), path("chunks/${chr_no}_all.txt"), emit: chunk_file
    tuple val(chr_no), path("tagged/${BIN}_${chr_no}.tagged.bcf*"), emit: bcf_tagged

    script:

    bcf_threads = (int) Math.ceil(task.cpus*2)
    """

    mkdir -p normalized
    mkdir -p tagged
    mkdir -p nogt
    mkdir -p chunks
    mkdir -p binary/${chr_no}_bin

    bcftools reheader -h ${reheader} ${vcf} --threads ${bcf_threads} | \
    bcftools norm -m -any -Ou --threads ${bcf_threads} | \
    bcftools view -m 2 -M 2 -v snps,indels --threads ${bcf_threads} -Ob -o normalized/${BIN}_${chr_no}.normalized.bcf
    bcftools index -f normalized/${BIN}_${chr_no}.normalized.bcf --threads ${bcf_threads}

    #adding allele number and allele count
    bcftools +fill-tags normalized/${BIN}_${chr_no}.normalized.bcf -Ob -o tagged/${BIN}_${chr_no}.tagged.bcf -- -t AN,AC
    bcftools index -f tagged/${BIN}_${chr_no}.tagged.bcf --threads ${bcf_threads}

    #create bcf with no GT infor for chunking
    bcftools view -G -Ob -o nogt/${BIN}_${chr_no}.nogt.bcf tagged/${BIN}_${chr_no}.tagged.bcf
    bcftools index -f nogt/${BIN}_${chr_no}.nogt.bcf --threads ${bcf_threads}

    #Chunking reference panel
    GLIMPSE2_chunk \
        --input nogt/${BIN}_${chr_no}.nogt.bcf \
        --region ${chr_no} \
        --output chunks/"${chr_no}_all.txt" \
        --map ${gmap} \
        --sequential

    """
}

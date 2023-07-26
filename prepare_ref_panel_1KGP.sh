mkdir original
mkdir tagged
mkdir nogt
mkdir chunks
mkdir binary
mkdir gmap
mkdir normalized

for i in {1..22}
do

#definition
in_vcf=CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz
gmap=chr${i}.b38.gmap.gz
norm_bcf=CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.normalized.bcf
tag_bcf=CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.tagged.bcf
nogt_bcf=CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.nogt.bcf
chunk=chr${i}_all.txt
BIN=CCDG_14151_B01_GRM_WGS_2020-08-05

#download reference panel from FTP
cd original
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${in_vcf}
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${in_vcf}.tbi
cd ../

#download genetic map
cd gmap
wget https://github.com/odelaneau/GLIMPSE/blob/master/maps/genetic_maps.b38/${gmap}
cd ../

#Normalization and filtering
bcftools norm -m -any original/${in_vcf} -Ou --threads 4 | \
bcftools view -m 2 -M 2 -V snps,indels --threads 4 -Ob -o normalized/${norm_bcf}
bcftools index -f normalized/${norm_bcf} --threads 4

#adding allele number and allele count
bcftools +fill-tags normalized/${norm_bcf} -Ob -o tagged/${tag_bcf} -- -t AN,AC
bcftools index -f tagged/${tag_bcf} --threads 4

#create bcf with no GT infor for chunking
bcftools view -G -Ob -o nogt/${nogt_bcf} tagged/${tag_bcf}
bcftools index -f nogt/${nogt_bcf} --threads 4

#Chunking reference panel
GLIMPSE2_chunk \
--input nogt/${nogt_bcf} \
--region chr${i} \
--output chunks/${chunk} \
--map gmap/${gmap} \
--sequential

#split reference panel into bin
mkdir binary/chr${i}_bin
while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)

GLIMPSE2_split_reference \
--reference tagged/${tag_bcf} \
--map gmap/${gmap} \
--input-region ${IRG} \
--output-region ${ORG} \
--output binary/chr${i}_bin/${BIN}
done < chunks/${chunk}

done

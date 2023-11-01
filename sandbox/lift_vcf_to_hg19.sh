#!/usr/bin/env bash

set -xeu

RESOURCE_DIR="/home/hw/Desktop/mincedTapestri/resources/"

vcf_file=$1

sed -i '2i##FORMAT=<ID=T,Number=1,Type=Integer,Description="?">' $vcf_file
sed -i '3i##FORMAT=<ID=E,Number=1,Type=Integer,Description="?">' $vcf_file
sed -i '4i##FORMAT=<ID=GO,Number=3,Type=Integer,Description="?">' $vcf_file
sed -i '5i##FORMAT=<ID=GN,Number=3,Type=Integer,Description="?">' $vcf_file
sed -i '6i##FILTER=<ID=BACKGROUND,Description="Background SNV">' $vcf_file

bcftools view -O z -o ${vcf_file}.gz ${vcf_file}
bcftools index ${vcf_file}.gz
# Remove partial chromosomes and pseudochromosomes
bcftools view \
	-R ${RESOURCE_DIR}/hg38.nuclearOnly.bed \
	-O z \
	-o ${vcf_file}.chromosomes.vcf.gz \
	${vcf_file}.gz

# Lift from hg38 to hg19
picard -Xmx12g LiftoverVcf \
	-I ${vcf_file}.chromosomes.vcf.gz \
	-O ${vcf_file}.chromosomes_hg19.vcf \
	-CHAIN ${RESOURCE_DIR}/hg38ToHg19.over.chain \
	-R ${RESOURCE_DIR}/hg19/hg19.fa.gz \
	-REJECT ${vcf_file}.chromosomes_rejected.vcf


bcftools view -O z -o ${vcf_file}.chromosomes_hg19.vcf.gz ${vcf_file}.chromosomes_hg19.vcf
bcftools index ${vcf_file}.chromosomes_hg19.vcf.gz
# Filter sites that overlap with target panel
bcftools view \
	-R ${RESOURCE_DIR}/Designer/3848-amplicon.bed \
	-O z \
	-o ${vcf_file}.final.vcf.gz \
	${vcf_file}.chromosomes_hg19.vcf.gz

rm ${vcf_file}.gz \
   ${vcf_file}.gz.csi \
   ${vcf_file}.chromosomes.vcf.gz \
   ${vcf_file}.chromosomes_hg19.vcf \
   ${vcf_file}.chromosomes_hg19.vcf.gz \
   ${vcf_file}.chromosomes_hg19.vcf.gz.csi \
   ${vcf_file}.chromosomes_hg19.vcf.idx \
   ${vcf_file}.chromosomes_rejected.vcf
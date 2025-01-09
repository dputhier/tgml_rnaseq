#!/bin/bash

mkdir -p input/fastq

cd input/fastq
rm -f *fastq.gz *fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/099/SRR11000299/SRR11000299.fastq.gz -O RipmOVA~mTEClo-int~SRR11000299.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/098/SRR11000298/SRR11000298.fastq.gz -O K14~mTEClo~int~SRR11000298.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/001/SRR11000301/SRR11000301.fastq.gz -O OTII~mTEClo~int~SRR11000301.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/095/SRR11000295/SRR11000295.fastq.gz -O CIITAKO~mTEClo~int~SRR11000295.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/094/SRR11000294/SRR11000294.fastq.gz -O mTEClo~int~rep~SRR11000294.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/096/SRR11000296/SRR11000296.fastq.gz -O CIITAKO~mTEClo~int~rep~SRR11000296.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/000/SRR11000300/SRR11000300.fastq.gz -O RipmOVA~mTEClo~int~rep~SRR11000300.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/002/SRR11000302/SRR11000302.fastq.gz -O OTII~mTEClo~int~rep~SRR11000302.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/093/SRR11000293/SRR11000293.fastq.gz -O mTEClo~int~SRR11000293.fq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/097/SRR11000297/SRR11000297.fastq.gz -O K14~mTEClo~int~SRR11000297.fq.gz

cd ..
mkdir -p  indexes
cd indexes
for i in /shared/bank/mus_musculus/mm10/star-2.7.5a/*; do ln -s $i; done

cd ..
mkdir -p gtf
cd gtf
cat  /shared/bank/mus_musculus/mm10/gff/Mus_musculus.GRCm38.97.gtf | grep -v "^#" | grep -P  "^[1-9MXY]" | sed 's/^/chr/' | sed 's/^chrMT/chrM/' > Mus_musculus.GRCm38.97.gtf

cd ..
mkdir -p  annotations
cd annotations
curl https://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/chromInfo.txt.gz | gunzip -c | grep -v "_"  > chromInfo_mm10.txt

cd ..
mkdir -p house_keeping
cd house_keeping
wget https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10.HouseKeepingGenes.bed.gz
gunzip -c mm10.HouseKeepingGenes.bed.gz | grep -v "chrX_GL456233_random" > mm10.HouseKeepingGenes.bed
rm -f mm10.HouseKeepingGenes.bed.gz

cd ..
mkdir -p fasta
cd fasta
ln -s /shared/bank/mus_musculus/mm10/fasta/mm10.fa

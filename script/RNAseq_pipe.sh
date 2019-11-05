#####RNA_seq#####
#1.Remove adaptor
bbduk.sh -Xmx16g in1=SRR2089613_1.fastq in2=SRR2089613_2.fastq out1=data_rm_adaptor_1.fq out2=data_rm_adaptor_2.fq ref=~/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
#2.Remove low quality reads
bbduk.sh -Xmx16g in1=data_rm_adaptor_1.fq in2=data_rm_adaptor_2.fq out1=data_clean_R1.fastq out2=data_clean_R2.fastq qtrim=r trimq=10 &
#3.mapping with hisat2
hisat2 -p 8 --dta-cufflinks -x ~/Annotation/hisat2_index/mm9_genome/genome -1 data_clean_R1.fastq -2 data_clean_R2.fastq -S data_RNAseq.sam
#4.samtools
samtools view -bS data_RNA.sam | samtools sort -@ 4 -o sort_data.bam
#5.cufflinks
cufflinks --no-update-check -p 8 -G ~/Annotation/ucsc_mm9_RNAseq.gtf -u --no-effective-length-correction --compatible-hits-norm -I 300000 -F 0.1 -j 0.15 -N -o data_out sort_data.bam
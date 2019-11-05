##!/bin/bash
#1.prepare TSS file
bedtools slop -i ~/Annotation/Chip_anno/Deeptools_input/UCSC/UCSC_mm9_TSS.bed -g ~/Annotation/Genome_anno/chrom_mm9.sizes -b 2500 > UCSC_mm9_slop_2.5kb.txt
#2.mapping
bowtie2 -p 8 -x ~/Annotation/Bowtie2_index/mm9/mm9 -1 data_H3K27ac_R1.fq.gz -2 data_H3K27ac_R2.fq.gz -S SAM/data_k27ac.sam
#3.samtools
samtools view -bS -F 4 data_k27ac.sam |samtools sort -@ 16 -o data_k27ac_sort.bam
#4.remove PCR duplicates
samtools rmdup data_k27ac_sort.bam data_k27ac_rmdup.bam
samtools index data_k27ac_rmdup.bam
#5.normalize data to 1x sequencing depth
bamCoverage -b data_k27ac_rmdup.bam -bs 1000 --normalizeUsing RPGC --effectiveGenomeSize 2304947926 --minMappingQuality 30 --ignoreForNormalization chrX chrY chrM -o data_H3K27ac_scale.bw
bigwigCompare -b1 data_H3K27ac_scale.bw -b2 data_input_scale.bw --pseudocount 1 --operation log2 -of bedgraph -o data_H3K27ac_ratio.bed
#6. Draw chip-seq figure
computeMatrix reference-point --referencePoint center -R 20kb_boundary.txt -b 600000 -a 600000 -S data_H3K27ac_scale.bw --skipZeros -bs 5000 -o k27ac_compute
plotProfile -m k27ac_compute -out k27ac.png
#7.calculate chip-seq signal in A/B compartment 
perl cal_mean_chipsig_region.pl data_A2B_compartment.txt data_H3K27ac_ratio.bed data_H3K27ac_a2b.txt
#8.get promoter(H3K4me3)
bedtools intersect -a data_k4me3.txt -b UCSC_mm9_slop_2.5kb.txt >promoter/promoters.bed
#9.get enhancer(H3K27ac-H3K27me3-H3K4me3 and > 2.5kb)
bedtools intersect -a data_k27ac.txt -b data_k27me3.txt -v | \
bedtools intersect -a - -b data_k4me3.txt -v | \
bedtools intersect -a - -b UCSC_mm9_slop_2.5kb.txt -v > enhancer/Enhancers.bed
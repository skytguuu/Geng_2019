#Hi-C part#
#bbmap clean data
bbduk.sh -Xmx32g in1=HiC1_R1.fastq.gz in2=HiC1_R2.fastq.gz out1=HiC1_bbmap2_R1.fastq.gz out2=HiC1_bbmap2_R2.fastq.gz ref=~/softwares/tools/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

#HiC-Pro#
HiC-Pro -i HiC_data -o HiC_result -c config_mm9.txt
#.hic file for Juicebox
bash ~/HiC-Pro/HiC-Pro_2.11.1/bin/utils/hicpro2juicebox.sh \
-i data.allValidPairs \
-g ~/HiC-Pro/HiC-Pro_2.11.1/annotation/chrom_mm9.sizes \
-j ~/juicer/juicer_tools.1.8.9_jcuda.0.8.jar \
-r ~/HiC-Pro/HiC-Pro_2.11.1/annotation/mm9_MboI.bed
#merge two replicates of Hi-C data
cat rep1.allValidPairs rep2.allValidPairs > merge_data.allValidPairs
#HiC-Pro for merge data
HiC-Pro -i merge_data/ -o merge_result/ -c config_mm9.txt -s build_contact_maps -s ice_norm
#merge .hic file
bash ~/HiC-Pro/HiC-Pro_2.11.1/bin/utils/hicpro2juicebox.sh \
-i merge_data.allValidPairs \
-g ~/HiC-Pro/HiC-Pro_2.11.1/annotation/chrom_mm9.sizes \
-j ~/juicer/juicer_tools.1.8.9_jcuda.0.8.jar \
-r ~/HiC-Pro/HiC-Pro_2.11.1/annotation/mm9_MboI.bed
#produce 40kb matrix
python ~/HiC-Pro-master/bin/utils/sparseToDense.py -b \
~/hic_results/matrix/merge/raw/40000/merge_40000_abs.bed \
~/hic_results/matrix/merge/iced/40000/merge_40000_iced.matrix \
-o 40kb --perchr

#DI calculate, TAD calling and chromatin loops prediction
#1.use 20kb matrix
python /share/home/Garen/Hicpro/HiC-Pro-master/bin/utils/sparseToDense.py -b \
~/hic_results/matrix/merge/raw/20000/merge_20000_abs.bed \
~/hic_results/matrix/merge/iced/20000/merge_20000_iced.matrix \
-o 20kb --perchr
#2.PSCYCHIC 
perl PSYCHIC_conf.pl
for i in /home/dell/mnt1/HiC_reference/development/HiTC/20kb/PSY/conf/zygote/*.conf
do 
python htad-chain.py $i 
done
cd Output
for i in *.domains.txt
do
name=$(basename $i)
chr=$(echo $name | cut -d _ -f1)
sed "s/^/$chr\t/g" $i > TAD/$i".txt"
done
#3.merge TAD file
cat chr*.txt > data_TAD.txt

#TAD signal
#1.change format 
for i in *.matrix
do 
paste row.txt $i | cat column.name - > $i.input 
done
#2.calculate insulation score#
perl ~/crane-nature-2015-master/scripts/matrix2insulation.pl -i merge_data_20kb_insulation.input -b 2000000 -ids 400000 -im mean -bmoe 3 -nt 0.1 -v

#statistic of cis interactions
perl static_intra-interaction_final_revise.pl merge_data.allValidPairs 20000 6 > merge_data_statistic.txt

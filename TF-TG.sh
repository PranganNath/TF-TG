#!/bin/bash 

run=B1

path=$(pwd)

mkdir -p $path/$run/ups/strand_+ve
mkdir -p $path/$run/ups/strand_-ve
mkdir -p $path/$run/ups/CRM


mkdir -p $path/$run/fimo/strand_+ve
mkdir -p $path/$run/fimo/strand_-ve
mkdir -p $path/$run/fimo/CRM
mkdir -p $path/$run/fimo/Merged


#Paths
pwms=/home/cluster/nath/Arabidopsis_1001/RESULTS/trial_05_10_23/pwms_all_motifs
list=/home/cluster/nath/Arabidopsis_1001/RESULTS/Benchmark_with_non_stress_responsive/SRGs/srg.txt
tf=/home/cluster/nath/Arabidopsis_1001/RESULTS/trial_05_10_23/TF_INFO.txt
ups=$path/$run/ups
pos=$path/$run/ups/strand_+ve
neg=$path/$run/ups/strand_-ve
crm=$path/$run/ups/CRM
fimo=$path/$run/fimo
fimo_pos=$path/$run/fimo/strand_+ve
fimo_neg=$path/$run/fimo/strand_-ve
fimo_crm=$path/$run/fimo/CRM
merged=$path/$run/fimo/Merged
ref_gen=/home/cluster/nath/Arabidopsis_1001/RESULTS/trial_05_10_23/Arabidopsis_thaliana.TAIR10.55.dna.toplevel.fa
matrix=/home/cluster/trainee/Arabidopsis_1001/trial_29-11-23/B1/alternate.R
matrix2=/home/cluster/trainee/Arabidopsis_1001/trial_29-11-23/B1/alternate2.R

#+ve stranded upstream genes
ups_coords_1=/home/cluster/nath/Arabidopsis_1001/DATA/ups/ups_pos_2798.txt

#-ve stranded upstream genes
ups_coords_2=/home/cluster/nath/Arabidopsis_1001/DATA/ups/ups_neg_2828.txt

#Genes containing CRMs
ups_coords_3=/home/cluster/nath/Arabidopsis_1001/DATA/ups/ups_crm_98.txt

details=/home/cluster/nath/Arabidopsis_1001/DATA/5626_genes/details_5626.txt

########################

:<< 'Upstream extraction using RSAT'

#STEP_1(For Docker): Create two directories for input files and results in working directory (wd)
mkdir -p rsat_data/genomes rsat_results
chmod -R a+w rsat_data/genomes rsat_results

#STEP_2(For Docker): Pull the docker RSAT container
docker pull biocontainers/rsat:20230828_cv1

#STEP_3(For Docker): Run the RSAT docker container
docker run --rm -v $wd/rsat_data:/packages/rsat/public_html/data/ -v $wd/rsat_results:/home/rsat_user/rsat_results -it biocontainers/rsat:20230828_cv1
##IGNORE APACHE ERROR

#STEP_4: Download the organism genomic data and files to be used 
download-organism -v 2 -org Arabidopsis_thaliana.TAIR10.55 -server https://rsat.eead.csic.es/plants

#STEP_5: Get the -1000 to -1 upstream regions of genes using RSAT
rsat retrieve-seq -org Arabidopsis_thaliana.TAIR10.55 -feattype gene -type upstream -format fasta -label id,name -from -1000 -to -1 -noorf -i Genes.txt -o RSAT.fa
##Genes.txt: Gene list (only IDs) RSAT.fa: Results fasta file

#STEP_6: Remove any warnings/unprocessed genes from the master set
grep -v "WARNING" RSAT.fa > ref_ups.fa

#STEP_7: Get the chromosome and positions of the upstream sequences from reference 
cat ref_ups.fa | grep '>' | awk -F ';' '{print $5}' | awk -F ':' '{print $3}' > chrom.txt
cat ref_ups.fa | grep '>' | awk -F ';' '{print $5}' | awk -F ':' '{print $4,$5,$6}' | tr ' ' '-' > positions.txt
paste chrom.txt positions.txt | tr '\t' ':' > ups.txt

#STEP_8: Categorize the genes into positive stranded and negative stranded
cat ups.txt | grep 'D' | sed 's/-D//g' > ups_pos.txt
cat ups.txt | grep 'R' | sed 's/-R//g' > ups_neg.txt

#STEP_9: Create Details file
paste <(cat ref_ups.fa | grep '>' | awk -F ';' '{print $1}' | cut -c 2- | tr '\t' ' ') <(cat ups.txt | sed 's/-D//g' | sed 's/-R//g') <(cat ups.txt) | tr '\t' ' ' > details.txt

#STEP_10: Plug the locations of the files generated to their respective variables:
#ups_coords_1=ups_pos.txt 
#ups_coords_2=ups_neg.txt 
#details=details.txt  

Upstream extraction using RSAT

#########################

echo "Indexing vcf files"
for m in $(ls *.vcf.gz); do tabix $m; done
echo "Indexing done"

########################
echo "Building upstream consensus from VCFs"

for n in $(ls *.vcf.gz)
do

for p in $(cat $ups_coords_1)
do
samtools faidx $ref_gen $p | bcftools consensus $n >> $n"_+ve".fa
done

for p in $(cat $ups_coords_2)
do
samtools faidx $ref_gen $p | bcftools consensus $n >> $n"_-ve".fa
done

for p in $(cat $ups_coords_3)
do
samtools faidx $ref_gen $p | bcftools consensus $n >> $n"_CRM".fa
done

mv *_+ve.fa $pos
mv *_-ve.fa $neg
mv *_CRM.fa $crm
#mv *.fa $ups
done
echo "Upstream consensus built"


########################
echo "Collecting all PWMs of TFs"

grep -f $list $tf > tflist.txt

while read -r line
do

	a=($line)
	motif="${a[0]}.txt"
	w=($(wc $pwms/$motif))
	if [[ $w -gt 1 ]]
	then

		echo "AC ${a[1]}" >> collected.txt
		echo "XX" >> collected.txt
		echo "ID ${a[2]}" >> collected.txt
		echo "XX" >> collected.txt
		cat $pwms/$motif >> collected.txt
		echo "XX" >> collected.txt
	fi
	
done < tflist.txt

echo "//" >> collected.txt

sed -i "s/Pos/PO/g" collected.txt
echo "All PWMs collected"
cp collected.txt $ups

########################
echo "PWM scanning of upstream regions"

cd $pos
for i in $(ls *.fa)
do

#Setting background with 1st order Markov Model

fasta-get-markov -m 1 -dna $i $i"_pos_background"

#Concatenating all motifs to one file (transfac2meme)

transfac2meme -bg $i"_pos_background" -use_acc $ups/collected.txt > $i"_pos_ref_seq.txt"

#Run pwm scanning 

fimo -o $fimo_pos/$i"_results" -bfile $i"_pos_background" --thresh 1.0E-5 $i"_pos_ref_seq.txt" $i

#Collecting only positive strand hits from FIMO
cat $fimo_pos/$i"_results"/fimo.tsv | awk '$5=="+" {print $0}' | awk '{print $1,$2,$3,$4,$5,$6,$7,$9}'  | awk '!x[$0]++' | tr '\t' ' ' > $fimo_pos/$i"_results"/$i"_hits.tsv"

#Generating edgelists
join <(sort -k2 $fimo_pos/$i"_results"/$i"_hits.tsv") <(sort -k3 $details) -1 2 -2 3 | awk '{print $2,$10,$9,$1,$11}' > $fimo_pos/$i"_results"/join1.txt
join <(sort -k1 $fimo_pos/$i"_results"/join1.txt) <(sort -k2 $details) -1 1 -2 2 | awk '{print $1,$6,$2,$3,$4}' > $fimo_pos/$i"_results"/edgelist.txt
done

cd $neg
for i in $(ls *.fa)
do
fasta-get-markov -m 1 -dna $i $i"_neg_background"
transfac2meme -bg $i"_neg_background" -use_acc $ups/collected.txt > $i"_neg_ref_seq.txt"
#sed -i 's/strands: + -/strands: -/g' $i"_pos_ref_seq.txt"
fimo -o $fimo_neg/$i"_results" -bfile $i"_neg_background" --thresh 1.0E-5 $i"_neg_ref_seq.txt" $i
cat $fimo_neg/$i"_results"/fimo.tsv | awk '$5=="-" {print $0}' | awk '{print $1,$2,$3,$4,$5,$6,$7,$9}' | awk '!x[$0]++' | tr '\t' ' ' > $fimo_neg/$i"_results"/$i"_hits.tsv"
join <(sort -k2 $fimo_neg/$i"_results"/$i"_hits.tsv") <(sort -k3 $details) -1 2 -2 3 | awk '{print $2,$10,$9,$1,$11}' > $fimo_neg/$i"_results"/join1.txt
join <(sort -k1 $fimo_neg/$i"_results"/join1.txt) <(sort -k2 $details) -1 1 -2 2 | awk '{print $1,$6,$2,$3,$4}' > $fimo_neg/$i"_results"/edgelist.txt
done

cd $crm
for i in $(ls *.fa)
do
fasta-get-markov -m 1 -dna $i $i"_crm_background"
transfac2meme -bg $i"_crm_background" -use_acc $ups/collected.txt > $i"_crm_ref_seq.txt"
fimo -o $fimo_crm/$i"_results" -bfile $i"_crm_background" --thresh 1.0E-5 $i"_crm_ref_seq.txt" $i
done

cd $fimo_pos
for i in $(for p in $(echo *_results); do echo $p | awk -F '_' '{print $1}'; done); do
cat $fimo_pos/$i*/edgelist.txt $fimo_neg/$i*/edgelist.txt >> $fimo/$i"_edgelist.txt"
done

#Creates unified network and plotting power law 
Rscript $matrix2 $fimo
 
cd $fimo
mv *.svg $merged
mv *_hubs.txt $merged
cat Merged_Network_Weighted_Edgelist_1.csv | awk -F ',' '$4>2 {print $0}' | awk '!x[$0]++' | tr ' ' ',' > $merged/Merged_Network_Weighted_Edgelist.csv
rm Merged_Network_Weighted_Edgelist_1.csv

#Finds most common element from the weighted edgelists (HCT counts)
most_common_element=$(cat $merged/Merged_Network_Weighted_Edgelist.csv | awk -F ',' '{print $1}' | sort | uniq -c | sort -k1,1nr | head -n1 | awk '{print $2}')
echo "Most,Common,Element:,$most_common_element" >> $merged/Merged_Network_Weighted_Edgelist.csv

#Plots power law for individual graphs
Rscript $matrix $fimo

cd $fimo
for i in $(ls *_matrix.csv)
do
mv -v $i "bak_"$i
echo -n "#", > $i
cat "bak_"$i >> $i
rm "bak_"$i
done

#Fixes filenames for random graphs
for i in $(ls *_random.svg); do mv $i $(echo $i | awk -F '_' '{print $1,$3,$5}' | tr ' ' '_'); done

#Fixes filenames for random superimposed graphs
for i in $(ls *si.svg); do mv $i $(echo $i | awk -F '_' '{print $1,$5}' | tr ' ' '_') ; done

#Fixes filenames for ecotype graphs
for i in $(ls *gr_svg); do mv $i $(echo $i | awk -F '_' '{print $1,$4}' | tr ' ' '.') ; done

#Fixes filenames for ecotype L-graphs
for i in $(ls *lg.svg); do mv $i $(echo $i | awk -F '_' '{print $1,$3}' | tr ' ' '_') ; done

#Fixes filenames for hub lists
for i in $(ls *_hubs.txt); do sed -i 's/"//g' $i; mv $i $(echo $i | awk -F '_' '{print $1,$3}' | tr ' ' '_') ; done

#Fixes filenames for weighted edgelists
for i in $(ls *weighted_edgelist.csv); do awk -F ',' '$7>2 {print $2,$4,$6,$7}' $i | awk '!x[$0]++' | tr ' ' ',' > $i"_edgelist.csv"; rm $i; done
for i in $(ls *_edgelist.csv); do mv $i $(echo $i | awk -F '_' '{print $1,$3,$5}' | tr ' ' '_'); done

#Finds most common element from the weighted edgelists (HCT counts)
for i in $(ls *weighted_edgelist.csv); do
most_common_element=$(cat $i | awk -F ',' '{print $1}' | sort | uniq -c | sort -k1,1nr | head -n1 | awk '{print $2}')
echo "Most,Common,Element:,$most_common_element" >> $i; done

#Fixes filenames for ecotype matrices
for i in $(ls *_matrix.csv); do mv $i $(echo $i | awk -F '_' '{print $1,$3}' | tr ' ' '_') ; done

echo "Pipeline run over!"
exit


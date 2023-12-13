#!/bin/bash  

fimo=/home/cluster/nath/Arabidopsis_1001/RESULTS/100_runs/A1/B1/B1/fimo/
fimo_pos=/home/cluster/nath/Arabidopsis_1001/RESULTS/100_runs/A1/B1/B1/fimo/strand_+ve
fimo_neg=/home/cluster/nath/Arabidopsis_1001/RESULTS/100_runs/A1/B1/B1/fimo/strand_-ve
matrix=/home/cluster/nath/Arabidopsis_1001/RESULTS/100_runs/alternate.R
matrix2=/home/cluster/nath/Arabidopsis_1001/RESULTS/100_runs/alternate2.R
merged=/home/cluster/nath/Arabidopsis_1001/RESULTS/100_runs/A1/B1/B1/fimo/Merged

#Setting ulimit for RAM utilization
ulimit -s unlimited

#mkdir -p $fimo/merged

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

#!/bin/bash
################################help#############################
## HPV whole-genome sequencing based on Illumina data
## $1 is the folder with the files to be processed
## $2 is the output folder

## OBS!!! Unfortunatly input and output folder need to be the same !!!!!

# Start up script:
# ./script.sh /path_to_input/ /path_to_output/ 


#########################################
########### Stats #######################
#########################################

cd $1
echo "filelist"
#ls -la *fastq.gz|cut -d" " -f10- > filelist.txt
ls *fastq.gz > filelist.txt

echo "count raw fastq files"
RAW="$1"/*fastq.gz
for r in $RAW
do

zcat $r|wc -l |awk '{print ($1/4)}'|sed 's/_L00[0-9]_R[0-9]_00[0-9].fastq.gz//g' > $r"wcl_raw.txt"
done

cat *"wcl_raw.txt" > wcl_raw.fq.txt
rm *"wcl_raw.txt"

echo "count trimmed files"
wc -l *001.fastq.gzadr_R*.fq|awk '{print ($1/4)}'|sed 's/_L00[0-9]_R[0-9]_00[0-9].fastq.gz//g' > wcl_adr_temp.txt
head -n -1 wcl_adr_temp.txt > wcl_adr.txt
rm wcl_adr_temp.txt

FILES="$1""/*001.fastq.gzadr_R*.fq"
echo "mean"
for f in $FILES
do

cat $f | awk '{if(NR%4==2) print length($1)}' | awk '{ total += $1 } END { print total/NR }' >$f"mean.txt"
done

cat *mean.txt > means.txt
rm *mean.txt

echo "print stats_reads_R1_R2.txt"
paste <(awk '{print $0}' filelist.txt) <(awk '{print $0}' wcl_raw.fq.txt) <(awk '{print $0}' wcl_adr.txt ) <(awk '{print $0}' means.txt ) > stats_temp.txt
grep "R1_001.fastq.gz" stats_temp.txt > R1_stats.txt
grep "R2_001.fastq.gz" stats_temp.txt > R2_stats.txt
paste <(awk '{print $0}' R1_stats.txt) <(awk '{print $0}' R2_stats.txt) > stats_reads_R1_R2_temp.txt
sort stats_reads_R1_R2_temp.txt | awk '{print $1,$2,$3,$4,$6,$7,$8}' stats_reads_R1_R2_temp.txt > stats_reads_R1_R2.txt

echo "print reads_aligned.txt"
grep "Assigned" *.fastq.gz.sam_HPV.txt.summary|sed 's/.sam_HPV.txt.summary//g'|sed 's/Assigned	//g'|sed 's/:/ /g'> reads_aligned_temp.txt
sort reads_aligned_temp.txt > reads_aligned.txt

echo "print reads_ontarget.txt"
grep "^gi_*" *.fastq.gz.sam_HPV.txt |sed 's/.sam_HPV.txt//g'|sed 's/:/\t/g' > reads_ontarget_temp.txt
awk '{a[$1]+=$8} END {for(i in a) print i" "a[i]}' reads_ontarget_temp.txt > reads_ontarget_unsorted.txt
sort reads_ontarget_unsorted.txt > reads_ontarget.txt

FILES="$1""/*.sam"
echo "insert size"
for f in $FILES
do

cat $f | awk '$9>0 && $9<1001 {print $9}' > $f"insertsize.txt"

done 

rename 's/S[0-9]*_L001_R1_001.fastq.gz.saminsertsize.txt/insertsize.txt/g' *

FILES="$1""/*insertsize.txt"
echo "mean insert size"
for f in $FILES
do
cat $f | awk '{ total += $1 } END { print total/NR }' > $f"mean.txt"

done 

awk '{print FILENAME, $1}' *mean.txt > insertmeans_temp.txt
mkdir "$1"/insertsize
mv *insertsize.txt "$1"/insertsize
rm *mean.txt

sed -i 's/_insertsize.txtmean.txt//g' insertmeans_temp.txt
sort insertmeans_temp.txt > insertmeans.txt

echo "print compiled_stats.txt"
paste <(awk '{print $1,$2,$3,$4,$5,$6,$7,$2,$3}' stats_reads_R1_R2.txt) <(awk '{print $2}' reads_aligned.txt) <(awk '{print $2}' reads_ontarget.txt) <(awk '{print $2}' insertmeans.txt) > compiled_stats_temp.txt
#paste <(awk '{print $1,$2,$3,$4,$5,$6,$7}' stats_reads_R1_R2.txt) <(awk '{print $3}' stats_reads_R1_R2.txt) <(awk '{print $2}' reads_aligned.txt) <(awk '{print $2}' reads_ontarget.txt) > compiled_stats_temp.txt
tr ' ' '\t' < compiled_stats_temp.txt  > compiled_stats.txt
sed -i 's/_S[0-9]*_L00[0-9]_R[0-9]_00[0-9].fastq.gz//g' compiled_stats.txt
 
echo "clean up"
rm filelist.txt
rm wcl_raw.fq.txt
rm wcl_adr.txt
rm means.txt
rm stats_reads_R1_R2_temp.txt 
rm stats_temp.txt
rm R1_stats.txt
rm R2_stats.txt
rm reads_aligned_temp.txt
rm reads_ontarget_temp.txt
rm reads_ontarget_unsorted.txt
rm compiled_stats_temp.txt
rm insertmeans_temp.txt
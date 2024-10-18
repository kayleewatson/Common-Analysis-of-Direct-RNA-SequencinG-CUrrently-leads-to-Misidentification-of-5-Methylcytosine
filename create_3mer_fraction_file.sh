#!/usr/bin/env bash

#this script requires 3 arguments:
## first argument is the path to the bed file of modified fractions (FILTERED_FILE)
## second argument is the path to the reference fasta
## third argument is the sample name

tmpfile=$(mktemp --suffix=.txt)
tmpfile2=$(mktemp --suffix=.txt)
flanking_tsv=$(mktemp --suffix=.tsv)
tmptsv=$(mktemp --suffix=.tsv)

awk '{print $1">"$2":"$3"<"$4"+"$5}' $1 > $tmpfile

# make sure flanking_nt_intervals.py script is in the same directory, or provide path
python flanking_nt_intervals.py $tmpfile > $tmpfile2

bedtools getfasta -fi $2 -bed $tmpfile2 -tab | awk -F ":" 'OFS="\t" {print $1,$2,$3}' | awk -F "-" 'OFS="\t" {print $1,$2,$3,$4}' > $tmptsv

join -j 2 -o 1.5,2.4 <(sort -k 1b,2 $tmpfile2) <(sort -k 1b,2 $tmptsv) | perl -pe 'tr/tT/uU/' | awk -v myvar=$3 '{print $1"\t"$2"\t"myvar}'

rm $tmpfile
rm $tmpfile2
rm $flanking_tsv
rm $tmptsv

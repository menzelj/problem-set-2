#! /usr/bin/env bash

data="../data-sets"

# Question 1
# Use BEDtools intersect to identify the size of the largest overlap
# between CTCF and H3K4me3 locations.

tfbs_bed=$data/bed/encode.tfbs.chr22.bed.gz
H3k4me3_bed=$data/bed/encode.h3k4me3.hela.chr22.bed.gz

zcat $tfbs_bed | awk '$4 == "CTCF"' > CTCF_peaks.bed

answer_1=$(bedtools intersect -a CTCF_peaks.bed -b $H3k4me3_bed -wo \
    | awk '{print $NF}' \
    | sort -nr | head -n1
)
                           
echo "answer-1: $answer_1" > answers.yml

# Question 2
# Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
# 19,000,500 on chr22 of hg19 genome build. Report the GC content as a
# fraction (e.g., 0.50).

gc_fasta=$data/fasta/hg19.chr22.fa

echo -e "chr22\t19000000\t19000500" > gc.bed

answer_2=$(bedtools nuc -fi $gc_fasta -bed gc.bed \
    | grep -v '^#' \
    | awk '{print $5}' )
echo "answer-2: $answer_2" >> answers.yml

# Question 3
# Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e.,
# interval) that has the largest mean signal in ctcf.hela.chr22.bg.gz.

signal=$data/bedtools/ctcf.hela.chr22.bg
peaks=$data/bed/peaks.chr22.bed
#| awk '$4 == "CTCF"'
#> peaks_CTCF.bed

answer_3=$(bedtools sort -i $peaks \
    | awk '$4 == "CTCF"' \
    | bedtools map -a stdin -b $signal -c 4 -o mean \
    | sort -k5n \
    | tail -n1 \
    | awk '{print $3 - $2}')
echo "answer-3: $answer_3" >> answers.yml

# Question 4
# Use BEDtools to identify the gene promoter (defined as 1000 bp
# upstream of
# a TSS) with the highest median signal in `ctcf.hela.chr22.bg.gz`.  Report
# the gene name (e.g., 'ABC123')

# mayby use tss.bed which already has sorted tss
# I may need to make a genome file with one interval of the whole chr22

echo -e "chr22\t52000000" > chr22_interval.bed

#tss.hg19.chr22.bed | awk 'BEGIN {OFS="\t"} ($6 == "+") \
#    {print $1, $2, $2 + 1, $4, $5, $6 }' > tss.hg19.chr22.plus.bed

answer_4=$(cat tss.hg19.chr22.bed | awk 'BEGIN {OFS="\t"} ($6 == "+") \ 
    {print $1, $2, $3, $4, $5, $6 }' \
    | bedtools flank -i stdin -g chr22_interval.bed -l 1000 -r 0 \
    > promoter_region.bed
cat tss.hg19.chr22.bed | awk 'BEGIN {OFS="\t"} ($6 == "-") \ 
    {print $1, $2, $3, $4, $5, $6 }' \
    | bedtools flank -i stdin -g chr22_interval.bed -l 0 -r 1000 \
    >> promoter_region.bed

bedtools sort -i promoter_region.bed > tmp_prom.bed
        mv tmp_prom.bed promoter_region.bed

bedtools map -a promoter_region.bed -b ctcf.hela.chr22.bg -c 4 -o median \
    | sort -k7n \
    | tail -n1 \
    | awk '{print $4}' )
echo "answer-4: $answer_4" >> answers.yml


## Question 5
# Use BEDtools to identify the longest interval on `chr22` that is not
# covered by `genes.hg19.bed.gz`. Report the interval like `chr1:100-500`.


answer_5=$(cat genes.hg19.bed \
    | bedtools complement -i stdin -g hg19.genome \
    | awk '{print $1, $2, $3, $3 -$2}' \
    | grep 'chr22' \
    | sort -k4n | tail -n1 \
    | awk '{print $1":"$2"-"$3}' )
echo "answer-5: $answer_5" >> answers.yml

# Question 6

answer_6=$( bedtools random -g hg19.genome | grep 'chr22'| head)
echo "answer-6: $answer_6" >> answers.yml

# Analysis of iso-seq reads: Fraction Spliced and Alternative Polyadenylation

# need parse_bam.py (modified from Schaerfen et al. 2022) and find_polyade.py

# map reads: iso-seq mode (max gap=intron 1kbp)
minimap2 -ax splice:hq -G 1k FI-OER-3-3.v2.6.fasta FI-OER-3-3_Ht_5_3.fasta.gz > iso-seq.aln.sam

# filter and sort reads (unmapped, not primary alignment, read fails platform/vendor quality checks, read is PCR or optical duplicate, supplementary alignment)
samtools view -F 3844 -@ 6 -Sb iso-seq.aln.sam | samtools sort -@ 6 -T xyz -o iso-seq.aln.filter.sort.bam
samtools index iso-seq.aln.filter.sort.bam 



## Calculate fraction spliced

# gff to bed
awk -F'\t' 'BEGIN {OFS = FS} {print $1,$4,$5,$9,$6,$7}' Hamiltosporidium_tvaerminnensis_FIOER33.v2.6.intron.gff > Hamiltosporidium_tvaerminnensis_FIOER33.v2.6.intron.strand.bed

bedtools bamtobed -bed12 -i iso-seq.aln.filter.sort.bam > iso-seq.aln.filter.sort.bed
bedtools intersect -wo -s -F 0.3 -a iso-seq.aln.filter.sort.bed -b Hamiltosporidium_tvaerminnensis_FIOER33.v2.6.intron.strand.bed | sort  > intron.filter.bed
# manually curated bed file with "real introns for which we have iso-seq coverage"
bedtools intersect -wo -s -F 0.3 -a iso-seq.aln.primary.bed -b Hamiltosporidium_tvaerminnensis_FIOER33.v2.6.intron.strand.good.bed | sort  > intron.good.spanning.filter.bed

paste <(awk '{print $14, $16}' intron.good.spanning.filter.bed | sort | uniq -c) <(awk '{if ($10>1) print $14, $16}' intron.good.spanning.filter.bed | sort | uniq -c)



## Alternative polyadenylation

# gff to bed
grep -v "#" FI-OER-3-3.v2.6/Hamiltosporidium_tvaerminnensis_FIOER33.gff3 | grep gene | awk -F'\t' 'BEGIN {OFS = FS} {print $1, $4, $5, $9, $6, $7}' > Hamiltosporidium_tvaerminnensis_FIOER33.v2.6.gene.bed

# get reads aligning to genes
samtools view -N <(bedtools intersect -wo -s -f 0.25 -abam iso-seq.aln.filter.sort.bam -b Hamiltosporidium_tvaerminnensis_FIOER33.v2.6.gene.bed -bed | cut -f4) -o iso-seq.aln.filter.sort.gene.bam iso-seq.aln.filter.sort.bam 

# get reads aligning to single genes
bedtools intersect -wo -s -f 0.25 -abam iso-seq.aln.filter.sort.bam -b Hamiltosporidium_tvaerminnensis_FIOER33.v2.6.gene.bed -bed > intersect.reads.genes.txt 
awk 'NR==FNR{seen[$4]++;next};seen[$4]==1' intersect.reads.genes.txt intersect.reads.genes.txt > intersect.unique.reads.genes.txt 

# get 3' ends
join -t $'\t' -1 4 -2 1 <(sort -k 4 intersect.unique.reads.genes.txt ) <(python parse_bam.py | sed -e 's/), (/\n/g' | sed "s/'//g" | sed 's/\[(//' | sed 's/)\]//' | sed 's/, /\t/g' | sort) | sort -k 16 > intersect.unique.reads.genes.stop.txt

python find_polyade.py

# martes
mkdir fasta_alignments

# extract read sequences from 40 nt upstream to 3' ends
while read -r line;
do
   strand=$(echo "$line" | cut -f2)
   scaff=$(echo "$line" | cut -f4)
   ID=$(echo "$line" | cut -f1 | sed -E 's/ID=([^;]*);.*/\1/')

   if [ $strand == "+" ]
   then
    p2=$(echo "$line" | cut -f5)
    p1=$((p2 - 40))
   else
       p1=$(echo "$line" | cut -f5)
    p2=$((p1 + 40))
   fi

   echo $strand
   echo $scaff
   echo $p1
   echo $p2
   
   java -jar ~/programs/jvarkit/dist/sam4weblogo.jar -r $scaff:$p1-$p2 iso-seq.aln.filter.sort.bam --naming %n | seqtk subseq - <(echo "$line" | cut -f6- | tr '\t' '\n' | sed "s/'//g") | sed -e 's/-//g' | sed -e 's/>>>>//g' | sed -e 's/>>>//g' | sed -e 's/>>//g' | sed -e 's/>[^m]//g' > tmp.fasta
   
   # reverse complement if necessary
   if [ $strand == "+" ]
   then
    cp tmp.fasta fasta_alignments/$ID.fasta
   else
       seqtk seq -r tmp.fasta | sed -e 's/>>>>//g' | sed -e 's/>>>//g' | sed -e 's/>>//g' | sed -e 's/>[^m]//g' > fasta_alignments/$ID.fasta
   fi

    rm tmp.fasta
   
done < peaks.highest.allinfo.500.tsv

# same for second highest peak
mkdir fasta_alignments_second

# extract read sequences from 40 nt upstream to 3' ends
while read -r line;
do
   strand=$(echo "$line" | cut -f2)
   scaff=$(echo "$line" | cut -f4)
   ID=$(echo "$line" | cut -f1 | sed -E 's/ID=([^;]*);.*/\1/')

   if [ $strand == "+" ]
   then
    p2=$(echo "$line" | cut -f5)
    p1=$((p2 - 40))
   else
       p1=$(echo "$line" | cut -f5)
    p2=$((p1 + 40))
   fi

   echo $strand
   echo $scaff
   echo $p1
   echo $p2
   
   java -jar ~/programs/jvarkit/dist/sam4weblogo.jar -r $scaff:$p1-$p2 iso-seq.aln.filter.sort.bam --naming %n | seqtk subseq - <(echo "$line" | cut -f6- | tr '\t' '\n' | sed "s/'//g") | sed -e 's/-//g' | sed -e 's/>>>>//g' | sed -e 's/>>>//g' | sed -e 's/>>//g' | sed -e 's/>[^m]//g' > tmp.fasta
   
   # reverse complement if necessary
   if [ $strand == "+" ]
   then
    cp tmp.fasta fasta_alignments_second/$ID.fasta
   else
       seqtk seq -r tmp.fasta | sed -e 's/>>>>//g' | sed -e 's/>>>//g' | sed -e 's/>>//g' | sed -e 's/>[^m]//g' > fasta_alignments_second/$ID.fasta
   fi

    rm tmp.fasta
   
done < peaks.second.allinfo.500.tsv

# same for third highest peak
mkdir fasta_alignments_third

# extract read sequences from 40 nt upstream to 3' ends
while read -r line;
do
   strand=$(echo "$line" | cut -f2)
   scaff=$(echo "$line" | cut -f4)
   ID=$(echo "$line" | cut -f1 | sed -E 's/ID=([^;]*);.*/\1/')

   if [ $strand == "+" ]
   then
    p2=$(echo "$line" | cut -f5)
    p1=$((p2 - 40))
   else
       p1=$(echo "$line" | cut -f5)
    p2=$((p1 + 40))
   fi

   echo $strand
   echo $scaff
   echo $p1
   echo $p2
   
   java -jar ~/programs/jvarkit/dist/sam4weblogo.jar -r $scaff:$p1-$p2 iso-seq.aln.filter.sort.bam --naming %n | seqtk subseq - <(echo "$line" | cut -f6- | tr '\t' '\n' | sed "s/'//g") | sed -e 's/-//g' | sed -e 's/>>>>//g' | sed -e 's/>>>//g' | sed -e 's/>>//g' | sed -e 's/>[^m]//g' > tmp.fasta
   
   # reverse complement if necessary
   if [ $strand == "+" ]
   then
    cp tmp.fasta fasta_alignments_third/$ID.fasta
   else
       seqtk seq -r tmp.fasta | sed -e 's/>>>>//g' | sed -e 's/>>>//g' | sed -e 's/>>//g' | sed -e 's/>[^m]//g' > fasta_alignments_third/$ID.fasta
   fi

    rm tmp.fasta
   
done < peaks.third.allinfo.500.tsv

mkdir xstreme

# run motif analysis
while read -r line;
do

   ID=$(echo "$line" | cut -f1 | sed -E 's/ID=([^;]*);.*/\1/')
    
   xstreme --minw 4 --rna --p fasta_alignments/$ID.fasta -o xstreme/$ID
done < peaks.highest.allinfo.500.tsv

cat fasta_alignments/*.fasta > all.sequences.fasta

xstreme --minw 4 --rna --p all.sequences.fasta -o xstreme.all

cat fasta_alignments_second/*.fasta > all.sequences_second.fasta

xstreme --minw 4 --rna --p all.sequences_second.fasta -o xstreme.all.second

cat fasta_alignments_third/*.fasta > all.sequences_third.fasta

xstreme --minw 4 --rna --p all.sequences_third.fasta -o xstreme.all.third



# extract ref sequences from 40 nt upstream to 3' ends and repeat motif analysis

mkdir ref_fasta

while read -r line;
do
   strand=$(echo "$line" | cut -f2)
   scaff=$(echo "$line" | cut -f4)
   ID=$(echo "$line" | cut -f1 | sed -E 's/ID=([^;]*);.*/\1/')

   if [ $strand == "+" ]
   then
    p2=$(echo "$line" | cut -f5)
    p1=$((p2 - 40))
   else
       p1=$(echo "$line" | cut -f5)
    p2=$((p1 + 40))
   fi

   echo $strand
   echo $scaff
   echo $p1
   echo $p2

   echo -e $scaff'\t'$p1'\t'$p2'\t'$ID'\t'"."'\t'$strand > ref_fasta/$ID.bed
   
   bedtools getfasta -s -fi FI-OER-3-3.v2.6.fasta -bed ref_fasta/$ID.bed >> ref_fasta/all.fasta
   
done < peaks.highest.allinfo.500.tsv

xstreme --minw 4 --rna --p ref_fasta/all.fasta -o xstreme.ref.first

# same for second highest peak
while read -r line;
do
   strand=$(echo "$line" | cut -f2)
   scaff=$(echo "$line" | cut -f4)
   ID=$(echo "$line" | cut -f1 | sed -E 's/ID=([^;]*);.*/\1/')

   if [ $strand == "+" ]
   then
    p2=$(echo "$line" | cut -f5)
    p1=$((p2 - 40))
   else
       p1=$(echo "$line" | cut -f5)
    p2=$((p1 + 40))
   fi

   echo $strand
   echo $scaff
   echo $p1
   echo $p2

   echo -e $scaff'\t'$p1'\t'$p2'\t'$ID'\t'"."'\t'$strand > ref_fasta/$ID.second.bed
   
   bedtools getfasta -s -fi FI-OER-3-3.v2.6.fasta -bed ref_fasta/$ID.second.bed >> ref_fasta/all_second.fasta
   
done < peaks.second.allinfo.500.tsv

xstreme --minw 4 --rna --p ref_fasta/all_second.fasta -o xstreme.ref.second

# same for third highest peak
while read -r line;
do
   strand=$(echo "$line" | cut -f2)
   scaff=$(echo "$line" | cut -f4)
   ID=$(echo "$line" | cut -f1 | sed -E 's/ID=([^;]*);.*/\1/')

   if [ $strand == "+" ]
   then
    p2=$(echo "$line" | cut -f5)
    p1=$((p2 - 40))
   else
       p1=$(echo "$line" | cut -f5)
    p2=$((p1 + 40))
   fi

   echo $strand
   echo $scaff
   echo $p1
   echo $p2

   echo -e $scaff'\t'$p1'\t'$p2'\t'$ID'\t'"."'\t'$strand > ref_fasta/$ID.third.bed
   
   bedtools getfasta -s -fi FI-OER-3-3.v2.6.fasta -bed ref_fasta/$ID.third.bed >> ref_fasta/all_third.fasta
   
done < peaks.third.allinfo.500.tsv

xstreme --minw 4 --rna --p ref_fasta/all_third.fasta -o xstreme.ref.third

# get position of motifs
while read -r line;
 do
 if [ `wc -l <(grep -obi 'A.TAAA' <<< "$line") | awk '{print $1}'` -ge "2" ]
 then
    grep -obi 'A.TAAA' <<< "$line" | grep -oE '[0-9]+' > tmp.pos
    awk -F"," -v c=1 -v t=20 '{a[NR]=$c}END{
        asort(a);d=a[NR]-t;d=d<0?-d:d;v = a[NR]
        for(i=NR-1;i>=1;i--){
                m=a[i]-t;m=m<0?-m:m
                if(m<d){
                    d=m;v=a[i]
                }
        }
        print v
    }' tmp.pos
 else
  grep -obi 'A.TAAA' <<< "$line" | grep -oE '[0-9]+'
 fi
done < ref_fasta/all.fasta > motif.pos.A.UAAA.txt

while read -r line;
 do
 if [ `wc -l <(grep -obi 'ATTAAA' <<< "$line") | awk '{print $1}'` -ge "2" ]
 then
    grep -obi 'ATTAAA' <<< "$line" | grep -oE '[0-9]+' > tmp.pos
    awk -F"," -v c=1 -v t=20 '{a[NR]=$c}END{
        asort(a);d=a[NR]-t;d=d<0?-d:d;v = a[NR]
        for(i=NR-1;i>=1;i--){
                m=a[i]-t;m=m<0?-m:m
                if(m<d){
                    d=m;v=a[i]
                }
        }
        print v
    }' tmp.pos
 else
  grep -obi 'ATTAAA' <<< "$line" | grep -oE '[0-9]+'
 fi
done < ref_fasta/all.fasta > motif.pos.AUUAAA.txt

while read -r line;
 do
 if [ `wc -l <(grep -obi 'AATAAA' <<< "$line") | awk '{print $1}'` -ge "2" ]
 then
    grep -obi 'AATAAA' <<< "$line" | grep -oE '[0-9]+' > tmp.pos
    awk -F"," -v c=1 -v t=20 '{a[NR]=$c}END{
        asort(a);d=a[NR]-t;d=d<0?-d:d;v = a[NR]
        for(i=NR-1;i>=1;i--){
                m=a[i]-t;m=m<0?-m:m
                if(m<d){
                    d=m;v=a[i]
                }
        }
        print v
    }' tmp.pos
 else
  grep -obi 'AATAAA' <<< "$line" | grep -oE '[0-9]+'
 fi
done < ref_fasta/all.fasta > motif.pos.AAUAAA.txt

while read -r line;
 do
 if [ `wc -l <(grep -obi 'TAAA' <<< "$line") | awk '{print $1}'` -ge "2" ]
 then
    grep -obi 'TAAA' <<< "$line" | grep -oE '[0-9]+' > tmp.pos
    awk -F"," -v c=1 -v t=18 '{a[NR]=$c}END{
        asort(a);d=a[NR]-t;d=d<0?-d:d;v = a[NR]
        for(i=NR-1;i>=1;i--){
                m=a[i]-t;m=m<0?-m:m
                if(m<d){
                    d=m;v=a[i]
                }
        }
        print v
    }' tmp.pos
 else
  grep -obi 'TAAA' <<< "$line" | grep -oE '[0-9]+'
 fi
done < ref_fasta/all.fasta > motif.pos.UAAA.txt

grep -c -i TAAA ref_fasta/all.fasta
grep -c -i A.TAAA ref_fasta/all.fasta
grep -c -i ATTAAA ref_fasta/all.fasta
grep -c -i AATAAA ref_fasta/all.fasta

grep -c -i TAAA ref_fasta/all_second.fasta
grep -c -i A.TAAA ref_fasta/all_second.fasta
grep -c -i ATTAAA ref_fasta/all_second.fasta
grep -c -i AATAAA ref_fasta/all_second.fasta

grep -c -i TAAA ref_fasta/all_third.fasta
grep -c -i A.TAAA ref_fasta/all_third.fasta
grep -c -i ATTAAA ref_fasta/all_third.fasta
grep -c -i AATAAA ref_fasta/all_third.fasta

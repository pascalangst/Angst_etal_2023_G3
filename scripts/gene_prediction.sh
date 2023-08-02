# PacBio ccs integration into gene prediction with BRAKER and TSEBRA
# see https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/long_reads/long_read_protocol.md

# following above tutorial:
stringtie2fa.py -g FI-OER-3-3.v2.6.fasta -f FI-OER-3-3.v2.6.collapsed.gff -o cupcake.v2.6.fa
gmst.pl --strand direct cupcake.v2.6.fa.mrna --output gmst.v2.6.out --format GFF
gmst2globalCoords.py -t FI-OER-3-3.v2.6.collapsed.gff -p gmst.v2.6.out -o gmst.v2.6.global.gtf -g FI-OER-3-3.v2.6.fasta

# convert to gff3
perl GeneMarkHMM_GTF_to_EVM_GFF3.pl gmst.v2.6.global.gtf > gmst.v2.6.global.gff3

# convert to oneline fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < FI-OER-3-3.v2.6.collapsed.rep.fa > FI-OER-3-3.v2.6.collapsed.rep.oneline.fasta
sed -i '/^$/d' FI-OER-3-3.v2.6.collapsed.rep.oneline.fasta

# predict genes
funannotate predict -i FI-OER-3-3.v2.6.fasta -o FI-OER-3-3.v2.6 -s "Hamiltosporidium tvaerminnensis" --isolate FIOER33 --busco_db microsporidia --busco_seed_species encephalitozoon_cuniculi_GB --cpus 10 --transcript_evidence FI-OER-3-3.v2.6.collapsed.rep.oneline.fasta --pasa_gff gmst.v2.6.global.gff3

# add UTRs and refine gene model predictions
funannotate update -i FI-OER-3-3.v2.6/ --pacbio_isoseq clustered.hq.fasta.gz --cpus 12

# Functional annotation
#No results from antiSMASH webserver
python iprscan-local.py -i FI-OER-3-3.v2.5/ -m local

funannotate annotate -i FI-OER-3-3.v2.6/ --cpus 12 --iprscan FI-OER-3-3.v2.6/annotate_misc/iprscan.xml --busco_db microsporidia

# update all genes annotated as "hypothetical proteins" with previous functional annotation
# get protein seqeunces
agat_sp_extract_sequences.pl -g Hamiltosporidium_tvaerminnensis_FIOER33.gff3 -f FI-OER-3-3.v2.6.fasta -p  -o infile_proteins.fa

# blast protein sequences to H. tvaerminnensis UniProt database
makeblastdb -in UP000292282.fasta -dbtype prot
blastp -query infile_proteins.fa  -db UP000292282.fasta -out protein.v1.blast -outfmt 6 -num_threads 16 # UniProt: Proteome ID UP000292282

# get all genes annotated as "hypothetical proteins" (n=2782) and their blast result
grep "=hypothetical protein;"  Hamiltosporidium_tvaerminnensis_FIOER33.gff3 | sed -E 's/.*ID=([^;]*);.*/\1/' > hypothetical.proteins.list
grep -f hypothetical.proteins.list protein.v1.blast > protein.v1.blast.hypothetical.proteins

# update gene and product name according to blast database
agat_sp_manage_functional_annotation.pl -f Hamiltosporidium_tvaerminnensis_FIOER33.gff3 -b protein.v1.blast.hypothetical.proteins --db  UP000292282.fasta -o Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3

# recover gene and product name from funannotate annotation where available
grep "product=hypothetical protein" Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3/Hamiltosporidium_tvaerminnensis_FIOER33.gff3  | sed -E 's/.*\t(ID[^;]*;).*/\1/' > overwritten.list
cat \
<( head -n 1 Hamiltosporidium_tvaerminnensis_FIOER33.gff3 ) \
<(grep $'\t'gene$'\t' Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3/Hamiltosporidium_tvaerminnensis_FIOER33.gff3 | grep -v -f <(sed -E 's/(ID=FUN_[0-9]*)-.*/\1/' overwritten.list)) \
<(grep $'\t'CDS$'\t' Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3/Hamiltosporidium_tvaerminnensis_FIOER33.gff3 ) \
<(grep $'\t'exon$'\t' Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3/Hamiltosporidium_tvaerminnensis_FIOER33.gff3 ) \
<(grep $'\t'five_prime_UTR$'\t' Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3/Hamiltosporidium_tvaerminnensis_FIOER33.gff3 ) \
<(grep $'\t'three_prime_UTR$'\t' Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3/Hamiltosporidium_tvaerminnensis_FIOER33.gff3 ) \
<(grep $'\t'mRNA$'\t' Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3/Hamiltosporidium_tvaerminnensis_FIOER33.gff3 | grep -v -f overwritten.list) \
<(grep $'\t'tRNA$'\t' Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3/Hamiltosporidium_tvaerminnensis_FIOER33.gff3 | grep -v -f overwritten.list) \
<(grep -f overwritten.list Hamiltosporidium_tvaerminnensis_FIOER33.gff3) \
<(grep $'\t'gene$'\t' Hamiltosporidium_tvaerminnensis_FIOER33.gff3 | grep -f <(sed -E 's/(ID=FUN_[0-9]*)-.*/\1/' overwritten.list)) \
> Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3/Hamiltosporidium_tvaerminnensis_FIOER33.high.quality.back.gff3

# sort gff
gff3sort.pl --precise Hamiltosporidium_tvaerminnensis_FIOER33.improved.gff3/Hamiltosporidium_tvaerminnensis_FIOER33.high.quality.back.gff3 > Hamiltosporidium_tvaerminnensis_FIOER33.lift-over.gff3

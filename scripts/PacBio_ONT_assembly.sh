## NextDenovo assembly
ls oxford_nanopore.fastq > input.fofn
nextDenovo run.cfg #genome_size = 20m # estimated genome size

# polishing 
ls FI-OER-3-3_R?_paired.fastq.gz > sgs.fofn
ls m54273_200909_122406.fasta.gz > lgs.fofn
nextPolish run_polish.cfg # Illumina and PacBio
medaka_consensus -i oxford_nanopore.fastq -d nd_polish/genome.nextpolish.fasta -o nd_polish_medaka -t 6 -m r941_min_high_g303 #ont errors

# polish nextDenovo assembly with pilon twice before merging:
minimap2 -ax map-pb -t 6 consensus.fasta pacbio.clr.canu.correctedReads.fasta.gz > consensus.pacbio.sam
samtools view -Sb consensus.pacbio.sam | samtools sort -T xyz -o consensus.pacbio.bam
rm consensus.pacbio.sam
samtools index -b consensus.pacbio.bam

bwa-mem2 index consensus.fasta
bwa-mem2 mem -t 12 consensus.fasta FI-OER-3-3_R1_paired.fastq.gz FI-OER-3-3_R2_paired.fastq.gz > consensus.sam
samtools view -@ 10 -Sb consensus.sam | samtools sort -@ 10 -T xyz -o consensus.bam
rm consensus.sam
samtools index -b consensus.bam

pilon --genome consensus.fasta --frags consensus.bam  --unpaired consensus.pacbio.bam --changes --output consensus.pilon -Xmx200G

# second round
minimap2 -ax map-pb -t 6 consensus.pilon.fasta pacbio.clr.canu.correctedReads.fasta.gz > consensus.pilon.pacbio.sam
samtools view -Sb consensus.pilon.pacbio.sam | samtools sort -T xyz -o consensus.pilon.pacbio.bam
rm consensus.pilon.pacbio.sam
samtools index -b consensus.pilon.pacbio.bam

bwa-mem2 index consensus.pilon.fasta
bwa-mem2 mem -t 12 consensus.pilon.fasta FI-OER-3-3_R1_paired.fastq.gz FI-OER-3-3_R2_paired.fastq.gz > consensus.pilon.sam
samtools view -@ 10 -Sb consensus.pilon.sam | samtools sort -@ 10 -T xyz -o consensus.pilon.bam
rm consensus.pilon.sam
samtools index -b consensus.pilon.bam

pilon --genome consensus.pilon.fasta --frags consensus.pilon.bam  --unpaired consensus.pilon.pacbio.bam --changes --output consensus.pilon2 -Xmx200G

# merge JFP's assembly (elongated PacBio assembly) with polished NextDenovo (nonopore) assembly using quickmerge
merge_wrapper.py Hami_FIOR33_2021_12_22.fasta consensus.pilon2.fasta -pre JFP_NextDenovo-np-medaka-pilon2

funannotate clean -i merged_JFP_NextDenovo-np-medaka-pilon2_manual3.fasta  -o merged_JFP_NextDenovo-np-medaka-pilon2_manual3.clean.fasta
funannotate sort -i merged_JFP_NextDenovo-np-medaka-pilon2_manual3.clean.fasta -o merged_JFP_NextDenovo-np-medaka-pilon2_manual3.sort.fa
funannotate mask -i merged_JFP_NextDenovo-np-medaka-pilon2_manual3.sort.fa -o merged_JFP_NextDenovo-np-medaka-pilon2_manual3.mask.fa
mv merged_JFP_NextDenovo-np-medaka-pilon2_manual3.mask.fa FI-OER-3-3.v2.6.fasta

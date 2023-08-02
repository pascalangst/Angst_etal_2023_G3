# load dependencies for karyoplot
library(karyoploteR)
library(seqinr)
library(rtracklayer)
library(ggplot2)
library(ggstance)

# generate custom genome
ff <- "FI-OER-3-3.v2.6.fasta"
fs <- read.fasta(file = ff)
Ht.genome <- toGRanges(data.frame(chr=gsub("scaffold_", "", getName(fs)), start=seq(1, 17), end=getLength(fs)))

# get repetitive regions
rep <- read.table("FI-OER-3-3.v2.6.fasta.out", quote="\"", comment.char="")
rep.range <- toGRanges(data.frame(chr=gsub("scaffold_", "", rep$V5), start=rep$V6, end=rep$V7))

# get telomeres
tel <- read.csv("telomeres_summary.csv")
tel$X....Contig <- gsub("scaffold_", "",tel$X....Contig)
tel$start <- ifelse(tel$X..of.rev.pattern..start. > 2, 6*tel$X..of.rev.pattern..start., 1)
tel$end <- ifelse(tel$X..of.pattern..end. > 2, as.numeric(gsub(",", "", tel$Length)) - 6*tel$X..of.pattern..end., as.numeric(gsub(",", "", tel$Length)))

tel.range <- toGRanges(data.frame(chr=c(tel$X....Contig, tel$X....Contig), start=c(rep(1,17),tel$end), end=c(tel$start,as.numeric(gsub(",", "", tel$Length)))))
tel.range <- tel.range[tel.range@ranges@width != 1]

# get buscos
busco <- read.delim("FI-OER-3-3.v2.6/full_table.tsv")
busco.range <- toGRanges(data.frame(chr=gsub("scaffold_", "",busco$Sequence), start=busco$Gene.Start, end=busco$Gene.End))

# make ideogram higher
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight = 100

# remove border and make black
trace(kpAddCytobands, edit = T)
#graphics::rect(xleft = xleft, xright = xright, ybottom = ybottom, 
#               ytop = ytop, col = "black", border = NA, ...)

# plot Karyotype
kp <- plotKaryotype(genome = Ht.genome, cytobands = rep.range, plot.params = pp)

# add telomeres, buscos, borders, and legend
kpRect(kp, data = tel.range, y0=0, y1=0.5, col="#FE6100", border="#FE6100", r0=0, r1=1)
kpRect(kp, data=busco.range, y0=0, y1=0.5, col="#648FFF", border="#648FFF")
kpRect(kp, data = Ht.genome, y0=-0.6, y1=-0.1, col=NA, border="gray", r0=0, r1=1)

legend(x = "bottomright", fill = c("#FE6100", "#648FFF", "black"), legend = c("Telomeres", "BUSCOs", "Repeats"))

# add axis
kpAddBaseNumbers(kp, tick.dist = 500000, minor.tick.dist = 250000,
                 add.units = T, cex=0.75, digits = 1)

# add all genes
gff.file <- "FI-OER-3-3.v2.6/Hamiltosporidium_tvaerminnensis_FIOER33.lift-over.gff3"
features <- import(gff.file)
features@seqnames <- gsub("scaffold_", "", features@seqnames)
genes <- features[features$type=="gene"]
kpAddCytobandsAsLine(kp)
kpPlotRegions(kp, data=genes)


# differentiated between long and short repeats to increase visibility
rep_long <- rep[rep$V7 - rep$V6 + 1 >= 1000,] 
rep_long.range <- toGRanges(data.frame(chr=gsub("scaffold_", "", rep_long$V5), start=rep_long$V6, end=rep_long$V7))
length(rep_long.range)

rep_short <- rep[rep$V7 - rep$V6 + 1 < 1000,] 
rep_short.range <- toGRanges(data.frame(chr=gsub("scaffold_", "", rep_short$V5), start=rep_short$V6, end=rep_short$V7))
length(rep_short.range)

trace(kpAddCytobands, edit = T)
#first_color <- "gray90"
#second_color <- "black"
#    for (i in seq_along(xleft)) {
#      if (i <= 17296) {
#        col <- first_color
#      }
#      else if (i <= 19246) {
#        col <- second_color
#      }
#      else {
#        col <- "green"
#      }
#      graphics::rect(xleft = xleft[i], xright = xright[i], 
#                     ybottom = ybottom[i], ytop = ytop[i], col = col, 
#                     border = NA, ...)
#    }

# plot Karyotype
kp <- plotKaryotype(genome = Ht.genome, cytobands = c(rep_short.range, rep_long.range), plot.params = pp)

# add telomeres, buscos, borders, and legend
kpRect(kp, data = tel.range, y0=0, y1=0.5, col="#FE6100", border="#FE6100", r0=0, r1=1)
kpRect(kp, data=busco.range, y0=0, y1=0.5, col="#648FFF", border="#648FFF")
kpRect(kp, data = Ht.genome, y0=-0.6, y1=-0.1, col=NA, border="gray", r0=0, r1=1)

legend(x = 0.78, y = 0.9, fill = c("#FE6100", "#648FFF", "black", "gray90"), legend = c("Telomeres", "BUSCOs", "Repeats \u2265 1 Kb", "Repeats < 1 Kb"))

# add axis
kpAddBaseNumbers(kp, tick.dist = 500000, minor.tick.dist = 250000,
                 add.units = T, cex=0.75, digits = 1)


# inset with repeats length profile
ggplot(rep, aes(x=V7 - V6 + 1)) + 
  geom_histogram(binwidth = 500, boundary = 0) + 
  geom_boxplot(aes(y=-500), width=750) + 
  theme_bw() + 
  ylab("Count") + xlab("Repeat length") + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15))

ggplot(rep, aes(x=V7 - V6 + 1)) + 
  geom_histogram() + 
  theme_bw() + 
  ylab("Count") + xlab("Repeat length") + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15)) + 
  scale_x_continuous(trans='log10', limits = c(1, max(rep$V7 - rep$V6 + 1)))

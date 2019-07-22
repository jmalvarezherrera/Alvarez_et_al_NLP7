library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
args = commandArgs(trailingOnly=TRUE)
bamfiles<-dir("./allbamfiles", recursive = TRUE, full.names = TRUE, pattern = "\\.sorted.bam$")
GFF_Loc<-file.path(args[1],"Araport11_GFF3_genes_transposons.201606.gff")
txdb_GFF<-makeTxDbFromGFF(GFF_Loc,format = 'gff')
ebg_GFF<-exonsBy(txdb_GFF,by="gene")
se<-summarizeOverlaps(features = ebg_GFF, reads = bamfiles,mode = "Union", ignore.strand = TRUE)
countscombined<-assay(se)
write.csv(countscombined, "counts_combined.csv")

#df_counts= read.csv(file='counts_combined.csv',header=TRUE)
#the_median <- apply(df_counts, 1, median)
#df_counts$Median <- the_median
#filtered_df_counts <- subset(df_counts, as.numeric(df_counts$Median)>10)
#filtered_df_counts$Median <- NULL
#write.csv(filtered_df_counts, "filtered_counts_combined.csv")
#cat('Total genes= ',nrow(df_counts))
#cat(', After removing low expressed genes = ',nrow(filtered_df_counts))


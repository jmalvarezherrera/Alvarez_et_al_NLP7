# DESeq2 analysis for TARGET experiment for TF2s
countscombined<-read.csv("counts_combined.csv",row.names = 1)

colnames(countscombined)<-c("EV_1","TGA4_1","CDF1_1","LBD37_1","LBD38_1",
                            "EV_2","TGA4_2","CDF1_2","LBD37_2","LBD38_2",
                            "EV_3","TGA4_3","CDF1_3","LBD37_3","LBD38_3")

groups<-as.factor(unlist(sapply(strsplit(colnames(countscombined),"_"),head,1)))
#the above is just a quick way to do the following if your names are formatted like I do them
#groups<-as.factor(c('"EV","TGA4","CDF1","LBD37","LBD38","EV","TGA4","CDF1","LBD37","LBD38","EV","TGA4","CDF1","LBD37","LBD38")))


aggregate<-as.matrix(countscombined)


###############################################################################################
# Filter low reads. Based on those with a median count that is less than 10 across all groups. 
###############################################################################################
grab<-function(a)
{output1<-sum(aggregate[a,1:c(ncol(aggregate))])}

maxes<-as.vector(NULL)
for (i in 1:length(aggregate[,1]))
{maxes[i] = grab(i)}


# Have a look at your distribution before filtering
x11()
hist(log(maxes,2),breaks=1000, xlim=c(1,15), ylim=c(0,300))


#determine the median for each set of replicates
medianCountByGroup = t(apply(aggregate, 1, tapply,
                             groups, median))
#calculte the max median from all samples for each gene.
maxMedian=apply(medianCountByGroup, 1, max)

#remove all genes where no group has a median # of counts of 10 or more
aggregate[maxMedian>=10,]->aggregate_filtered

# Get the sum of counts for each row after filtering
grab<-function(a)
{output1<-sum(aggregate_filtered[a,1:ncol(aggregate)])}

maxes_filtered<-as.vector(NULL)
for (i in 1:length(aggregate_filtered[,1]))
{maxes_filtered[i] = grab(i)}


# Have a look at your distribution after filtering
x11()
hist(log(maxes_filtered,2),breaks=1000, xlim=c(1,15), ylim=c(0,100))

#Remove rows that are not genes (e.g. ambiguous read counts, RFP/GFP, Vector).
aggregate_filtered<-aggregate_filtered[!rownames(aggregate_filtered) %in% c("5SrRNA","18SrRNA","25S/18S","5.8S/25S/18rRNA","__ambiguous","__alignment_not_unique","__no_feature","__not_aligned","__too_low_aQual","Vector_start","Vector_intermediate","Vector_end","RFP/mCherry","GR","GFP"),]

###############################################################################################
# do the EDA normalisation, only used for making MDS plot and clustering dendogram
# EDASeq package is not installed on mercer
###############################################################################################
library(EDASeq)
#Normalization TMM using EdgeR
library(edgeR)
dge<-DGEList(aggregate_filtered)
dge<-calcNormFactors(dge)
fullmat<-dge$counts

#Make MDS Plot (not making the pdf for some reason when run in the script, but works when run separataly)

expression.dist<-dist(t(fullmat[,1:ncol(fullmat)]))
expression.mds<-cmdscale(expression.dist,eig=TRUE, k=2)

x.mds<-expression.mds$points[,1]
y.mds<-expression.mds$points[,2]
temp<-as.data.frame(cbind(x.mds,y.mds))
row.names(temp)<-colnames(aggregate)

library(ggplot2)
x11()
ggplot(temp,aes(x=x.mds,y=y.mds,color=groups,label=rownames(temp)))+
  geom_point()+
  geom_text(hjust=0,nudge_x = 0.5) +
  theme(panel.background = element_rect(fill ="white",color="black"))
ggsave('MDSPlot.pdf')

pdf("Clustering.pdf", width=8, height=8)
normalized.clust = hclust(as.dist(1-cor(fullmat)), "ave")
print(plot(normalized.clust,xlab="Heirarchical Clustering using Pearson Correlation"))
dev.off()

###############################################################################################
#Differential Expression using DESeq2
###############################################################################################

#reorder the matrix, get rid of extra columns
EVcol<-grep("EV",colnames(fullmat))

TGA4col<-grep("TGA4",colnames(fullmat))
CDF1col<-grep("CDF1",colnames(fullmat))
LBD37col<-grep("LBD37",colnames(fullmat))
LBD38col<-grep("LBD38",colnames(fullmat))


newmat<-aggregate_filtered[,c(EVcol,TGA4col,CDF1col,LBD37col,LBD38col)]
TFgroups<-groups[c(EVcol,TGA4col,CDF1col,LBD37col,LBD38col)]
TFgroups<-relevel(TFgroups,ref="EV")

TFmeta<-data.frame(Treatment=TFgroups,row.names = colnames(newmat))

library(DESeq2)

DEmat<-DESeq(DESeqDataSetFromMatrix(countData = newmat,colData = TFmeta,design = ~Treatment))
resTGA4<-results(DEmat,contrast = c("Treatment","TGA4","EV"))
resCDF1<-results(DEmat,contrast = c("Treatment","CDF1","EV"))
resLBD37<-results(DEmat,contrast = c("Treatment","LBD37","EV"))
resLBD38<-results(DEmat,contrast = c("Treatment","LBD38","EV"))

#subset the results at FDR 0.1 and write to file.
write.csv(subset(resTGA4,padj<0.1),"TGA4-EV.csv")
write.csv(subset(resCDF1,padj<0.1),"CDF1-EV.csv")
write.csv(subset(resLBD37,padj<0.1),"LBD37-EV.csv")
write.csv(subset(resLBD38,padj<0.1),"LBD38-EV.csv")


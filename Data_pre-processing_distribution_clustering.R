##########################################################################################################################
### PART 1) DATA PRE-PROCESSING
##########################################################################################################################

#set the current folder as the working directory 
setwd("~/");

#data input
MyData <- read.csv(file="275genes_vs_samples.csv", header=TRUE, sep=",");

#remove the first column with sample-which is just the sample names repeating again
MyData2 <- MyData[,-c(1)];

#make the first column (samples names) as row names
MyData2_2 <- data.frame(MyData2[,-1], row.names=MyData2[,1])

#remove the rows that do NOT have full data (in this case,from 736 rows,only 434 will remain)
MyData3 <- na.omit(MyData2_2);
 
#transpose the matrix for setting genes as rows (for a further clustering by the expression profile)
MyData4 <- t(MyData3);

#########################################################################################################################
##### PART 2) CHECK THE DATA DISTRIBUTION
#########################################################################################################################

##plot the data to check the distribution (as doing correlation to very skewed data can return spurious results)
d <- density(MyData4)

plot (d, main="Figure 1. Distribution of data across TCGA colorectal cancer dataset", #add a title with main
      yaxs="i", #remove the small space between X-axis and the distribution area with y-axis="i"
      xlab="Samples;  N=119,350 elements",#add x-axis which are the samples & no of samples (119350elements=434reads for each one of 275genes)
      sub=expression(paste(Log[2], "(x+1), where x is RSEM normalised")), #add a subtitle with sub function, explaining which units are used for the samples
      mgp=c(1.9,.5,0) )#adjust axis label position with respect to the tick lables so they are at the right distance with function mgp

polygon(d, col="blue", border = "black"); #color the distribution area under the curve for a better visualisation 
#y-axis (i.e. density) is calculated by the probability density function

#########################################################################################################################
### PART 2) CLUSTERING BY GENES & BY SAMPLES IN A HEATMAP
#########################################################################################################################

##cluster genes by expression profile-looking for enriched/distinct groups by calculating the absolute correlation distance & then organise them on a heatmap for visualisation
##the absolute correlation used as it may be possible for some genes to be over-expressed, while other under-expressed - which may show results for canonical vs non-canonical
#1) pairwise correlation matrix between values on row and on column; method is Pearson correlation
#2) convert correlation matrix as a distance matrix by 1-|r| (where r=Pearson correlation coefficient)
#3) clustering is done with the results of the distance matrix

#1. row-clustering (hr) and column-wise clustering (hc)
#2. converting correlation matrix into a distance matrix 
#3. hierarchical clustering with dendrogram
hr <- hclust(as.dist(1-cor(t(MyData4), method="pearson")), method="complete")
#as scale function will be used (which works by column), the matrix has to be transposed as the genes are represented as rows
hc <- hclust(as.dist(1-cor(MyData4, method="pearson")), method="complete") 
##Cutting the dendrogram with function cutree() for estabilishing a number of clusters by gene
#the dendrogram of genes is cut at a specific height i.e. h() -by dividing by 1.5 each distance tree
#the clusters are divided by a color column vector in rainbow colours, so the differences between clusters can be spotted easily
genes_dendrogram_cut <- cutree(hr, h=max(hr$height)/1.5); 
mycolhc <- rainbow(length(unique(genes_dendrogram_cut)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(genes_dendrogram_cut)] 
#cut the dendrogram for samples and divide clusters in rainbow colours
samples_dendrogram_cut <- cutree(hc, h=max(hc$height)/1.5); 
myrowhc <- rainbow(length(unique(samples_dendrogram_cut)), start=0.1, end=0.9); myrowhc <- myrowhc[as.vector(samples_dendrogram_cut)]

#Plot heatmap 
library(gplots)#gplots package is loaded as it is required for using heatmap.2() function
library(RColorBrewer)#RColorBrewer package is loaded as it helps in generating a color palette that can help in visualising the clusters in the data


colour_palette  <- brewer.pal(4, "RdYlBu")
##set the colour palette as diverging from Red to Yellow to Green with 4 types of colours, as the main interest is in the middle two colours, 
#which represent most of the data (between 4 and 12)
##create a heatmap with the hierarchically clustered rows and columns from above and set the colour from the colour_palette from above
#the rows will be set as the scale
#with xlab and ylab, axes labels will be added
#since the size of the graph does not allow reading the name of each gene and each sample, these can be removed with labCol and labRow
#margins between axes labels and the actual plot are adjusted with margins() function
#copy the plot and save it as a .png format with dev.copy() & dev.off(); the plot will be found in the current working directory
#the heatmap is based on MyData4
#the dendrograms for rows and columns are shown on the figure, after they were previously cut at a specific distance
#no density and trace information is added on the graph for the easiness of interpretation
#by scale="row" argument, each value is converted into a Z-score for the row
dev.copy(png, 'Figure 2. Correlation heatmap of genes versus samples across colorectal cancer data from TCGA project.png')
heatmap.2 (MyData4, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=colour_palette, scale="row", 
           density.info="none", trace="none", RowSideColors=mycolhc, ColSideColors=myrowhc, 
           xlab="Samples", ylab="Genes", labRow=FALSE, labCol=FALSE, 
           key=TRUE,
           margins=c(1.5,1.5),
           main="Figure 2. Correlation heatmap of genes versus samples 
across colorectal cancer data from TCGA project")
dev.off()

#return the genes (rows) after the hierarchical clustering, so they can be assessed separately 
genes_after_clustering <- (MyData4[rev(hr$labels[hr$order]), hc$labels[hc$order]])[,0]
genes_after_clustering
#write the genes in the "hierarchical order" to a text file; the file will be found in the current working directory
write.table(genes_after_clustering, file = "genes_after_clustering.txt", quote=FALSE, sep="")

#look at the cluster number assigned to each gene
cluster_dendrogram_cut

#organise the genes by clusters
genes_by_clusters <- cluster_dendrogram_cut[hr$order]
genes_by_clusters

#compute the no of genes in each cluster
nrow(genes_after_clustering[1:14,])
nrow(genes_after_clustering[15:20,])
nrow(genes_after_clustering[21:41,])
nrow(genes_after_clustering[42:51,])
nrow(genes_after_clustering[52:85,])
nrow(genes_after_clustering[86:132,])
nrow(genes_after_clustering[133:148,])
nrow(genes_after_clustering[149:162,])
nrow(genes_after_clustering[163:173,])
nrow(genes_after_clustering[174:193,])
nrow(genes_after_clustering[194:199,])
nrow(genes_after_clustering[200:216,])
nrow(genes_after_clustering[217:225,])
nrow(genes_after_clustering[226:245,])
nrow(genes_after_clustering[246:259,])
nrow(genes_after_clustering[260:275,])


#extract the first x genes, so we can asses gene names to different clusters for further analysis
#save the genes from each cluster (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16) to a text file with the no of the cluster and no of genes inside that specific cluster
write.table(genes_after_clustering[1:14,], file="cluster1_14genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[15:20,], file="cluster2_6genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[21:41,], file="cluster3_21genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[42:51,], file="cluster4_10genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[52:85,], file="cluster5_34genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[86:132,], file="cluster6_47genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[133:148,], file="cluster7_16genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[149:162,], file="cluster8_14genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[163:173,], file="cluster9_11genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[174:193,], file="cluster10_20genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[194:199,], file="cluster11_6genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[200:216,], file="cluster12_17genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[217:225,], file="cluster13_9genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[226:245,], file="cluster14_20genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[246:259,], file="cluster15_14genes.txt", quote=FALSE, sep="")
write.table(genes_after_clustering[260:275,], file="cluster16_16genes.txt", quote=FALSE, sep="")

#extract the names of the column clusters (i.e. clusters of samples) to check if there is a dependence relationship
samples_after_clustering <- (MyData4[rev(hr$labels[hr$order]), hc$labels[hc$order]])[0,]

#write the genes in the "hierarchical order" to a text file
write.table(samples_after_clustering, file = "samples_after_clustering.txt", quote=FALSE, sep="")

#save the clusters by samples (A,B,C)
#samples_after_clustering will show samples as columns, so t() i.e. transpose function will convert them to rows
write.table(t(samples_after_clustering[,1:244]), file="clusterA_244samples.txt", quote=FALSE, sep="")
write.table(t(samples_after_clustering[,245:381]), file="clusterB_137samples.txt", quote=FALSE, sep="")
write.table(t(samples_after_clustering[,382:434]), file="clusterC_53 samples.txt", quote=FALSE, sep="")

##(if needed) Extracting a specific gene (row) with the samples
mean(MyData4 ["FZD9",])
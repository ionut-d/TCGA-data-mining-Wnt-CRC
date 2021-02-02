##########################################################################################################################
### PART 3) KRUSKAL-WALLIS TEST (for clustering by sample)
##########################################################################################################################

#read the full clinical data file as a table (because it is .tsv format)
clinical_data <- read.table(file="clinical.data.tsv", sep="\t", header=TRUE)

#eliminate the first column (which is just repeating itself)
data_1 <- clinical_data[,-c(1)]

####CLUSTER_A OF SAMPLES
#upload clusterA as a matrix
cluster_A <- read.table(file="clusterA_244samples.txt", header=FALSE, sep="")
#select only clinical data for rows (i.e. samples) from cluster A
data_3 <- merge(data_1, cluster_A, by=1:1)
#remove the gene expression columns and keep only the clinical information columns
data_5 <- data_3[,-c(2:276)]
#look at the classes of variables
str(data_5)
summary(data_5$sample_type)
#resolve the problem of missing data (replace NULL values with NA)
data_5[data_5==""] <- "No value"
#see summary of data_5 dataframe (i.e. how many in each category)
summary(data_5$tissue_source_site)
print(chisq.test(data_5$pathologic_stage))

####CLUSTER_B OF SAMPLES
#upload clusterB as a matrix
cluster_B <- read.table(file="clusterB_137samples.txt", header=FALSE, sep="")
#select only clinical data for rows (i.e. samples) from cluster A
data_6 <- merge(data_1, cluster_B, by=1:1)
#remove the gene expression columns and keep only the clinical information columns
data_7 <- data_6[, -c(2:276)]
#look at the classes of variables
str(data_7)
summary(data_7$sample_type)
#resolve the problem of missing data (replace NULL values with NA)
data_7[data_7==""] <- "No value"
#see summary of data_7 dataframe (i.e. how many in each category)
summary(data_7$tissue_source_site)

####CLUSTER_C OF SAMPLES
#upload clusterC as a matrix
cluster_C <- read.table(file="clusterC_53 samples.txt", header=FALSE, sep="")
#select only clinical data for rows (i.e. samples) from cluster A
data_8 <- merge(data_1, cluster_B, by=1:1)
#remove the gene expression columns and keep only the clinical information columns
data_9 <- data_8[, -c(2:276)]
#look at the classes of variables
str(data_9)
summary(data_9$sample_type)
#resolve the problem of missing data (replace NULL values with NA)
data_9[data_9==""] <- "No value"
#see summary of data_7 dataframe (i.e. how many in each category)
summary(data_9$tissue_source_site)

#test for sample type (result: null hypothesis)
x <- summary(data_5$sample_type)
y <- summary(data_7$sample_type)
z <- summary(data_9$sample_type)
kruskal.test(list(x, y, z))

#test for pathologic stage (result: null hypothesis)
x <- summary(data_5$pathologic_stage)
y <- summary(data_7$pathologic_stage)
z <- summary(data_9$pathologic_stage)
kruskal.test(list(x, y, z))

#test for tissue source site(result: p=0.0535); do a post-hoc test to see which group is different
x <- summary(data_5$tissue_source_site)
y <- summary(data_7$tissue_source_site)
z <- summary(data_9$tissue_source_site)
kruskal.test(list(x, y, z))

#do a post-hoc test to see which group is different (Dunn's test) and apply correction of p-value Benjamini-Hochberg
library(dunn.test)
p.adjustment.methods=c("bh")
dunn.test(list(x,y,z), method=p.adjustment.methods, kw=TRUE, label=TRUE, wrap=FALSE, table=TRUE, list=FALSE, alpha=0.05, altp=FALSE)

#result=clusterA different to clusterB (maybe colon cancer vs rectal cancer?)
#do a hypergeometric test to see what were the chances to actually get a colon cancer (successes) from the total no of samples (trials)
#hypergeometric is without replacement
#probability of getting colon cancer in the 244 samples from cluster_A
#lower.tail=FALSE is used for under-representation as we want to know the probability to right, so at least or equal with X (i.e. P>/= to x, what we are interested); x or more successes
#lower.tail=FALSE give area to the right; same as 1-phyper(..,lower.tail=TRUE,)
phyper(186-1, 286, 95, 244, lower.tail=TRUE, log.p=FALSE) #result=0.28
#probability of getting rectal cancer in the 244 samples from cluster_A
phyper(58-1, 95, 286, 244, lower.tail=TRUE, log.p=FALSE) #result=0.80

#result=clusterA different to clusterB (maybe rectal cancer vs colon cancer?)
#do a hypergeometric test to see what were the chances to actually get a rectal cancer (successes) from the total no of samples (trials)
#probability of getting colon cancer in the 137 samples from cluster_B; result=0.28
phyper(100-1, 286, 95, 137, lower.tail=TRUE, log.p=FALSE) #result=0.80
#probability of getting rectal cancer in the 137 sammples from cluster_B; result=0.80
phyper(37-1, 95, 286, 137, lower.tail=TRUE, log.p=FALSE) #result=0.28


#test for histlogical type (result: null hypothesis)
x <- summary(data_5$histological_type)
y <- summary(data_7$histological_type)
z <- summary(data_9$histological_type)
kruskal.test(list(x, y, z))

#test for anatomical neoplasm subdivision (result: null hypothesis)
x <- summary(data_5$anatomic_neoplasm_subdivision)
y <- summary(data_7$anatomic_neoplasm_subdivision)
z <- summary(data_9$anatomic_neoplasm_subdivision)
kruskal.test(list(x, y, z))

#test for tumor tissue site (result: null hypothesis)
x <- summary(data_5$tumor_tissue_site)
y <- summary(data_7$tumor_tissue_site)
z <- summary(data_9$tumor_tissue_site)
kruskal.test(list(x, y, z))

#test for follow-up for treatment succes (result: null hypothesis)
x <- summary(data_5$followup_treatment_success)
y <- summary(data_7$followup_treatment_success)
z <- summary(data_9$followup_treatment_success)
kruskal.test(list(x, y, z))

#test for colon polyps present (result: null hypothesis)
x <- summary(data_5$colon_polyps_present)
y <- summary(data_7$colon_polyps_present)
z <- summary(data_9$colon_polyps_present)
kruskal.test(list(x, y, z))

#test for a history of colon polyps (result: null hypothesis)
x <- summary(data_5$history_of_colon_polyps)
y <- summary(data_7$history_of_colon_polyps)
z <- summary(data_9$history_of_colon_polyps)
kruskal.test(list(x, y, z))

#test for history of neoadjuvant treatment (result: chi-squared value of 0, because none of the groups have history of neoadjuvant treatment)
x <- summary(data_5$history_of_neoadjuvant_treatment)
y <- summary(data_7$history_of_neoadjuvant_treatment)
z <- summary(data_9$history_of_neoadjuvant_treatment)
kruskal.test(list(x, y, z))

#test from KRAS mutations (result: null hypothesis)
x <- summary(data_5$kras_mutation_found)
y <- summary(data_7$kras_mutation_found)
z <- summary(data_9$kras_mutation_found)
kruskal.test(list(x, y, z))

#test for non-nodal tumour deposits (result: null hypothesis)
x <- summary(data_5$non_nodal_tumor_deposits)
y <- summary(data_7$non_nodal_tumor_deposits)
z <- summary(data_9$non_nodal_tumor_deposits)
kruskal.test(list(x, y, z))

#test for number of abnormal loci (result: null hypothesis)
x <- summary(data_5$number_of_abnormal_loci)
y <- summary(data_7$number_of_abnormal_loci)
z <- summary(data_9$number_of_abnormal_loci)
kruskal.test(list(x, y, z))

#test for sample type ID (result:chi-squared very small=0.049927, maybe reject null hypothesis?)
x <- summary(data_5$sample_type_id)
y <- summary(data_7$sample_type_id)
z <- summary(data_9$sample_type_id)
kruskal.test(list(x, y, z))

#do a post-hoc test to see which group is different (Dunn's test) and apply correction of p-value Benjamini-Hochberg
library(dunn.test)
p.adjustment.methods=c("bh")
dunn.test(list(x,y,z), method=p.adjustment.methods, kw=TRUE, label=TRUE, wrap=FALSE, table=TRUE, list=FALSE, alpha=0.05, altp=FALSE)

#test for synchronous colon cancer presence (result: null hypothesis)
x <- summary(data_5$synchronous_colon_cancer_present)
y <- summary(data_7$synchronous_colon_cancer_present)
z <- summary(data_9$synchronous_colon_cancer_present)
kruskal.test(list(x, y, z))

#test for venous invasion (result: null hypothesis)
x <- summary(data_5$venous_invasion)
y <- summary(data_7$venous_invasion)
z <- summary(data_9$venous_invasion)
kruskal.test(list(x, y, z))

#test for new tumour event after initial treatment (result: null hypothesis)
x <- summary(data_5$new_tumor_event_after_initial_treatment)
y <- summary(data_7$new_tumor_event_after_initial_treatment)
z <- summary(data_9$new_tumor_event_after_initial_treatment)
kruskal.test(list(x, y, z))

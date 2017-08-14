#Determine whether or not case status is associated with average values across different features (e.g. are CpG islands differentially methylated)

library('CpGassoc')
library(minfi)


#get number of rows in file (reads file faster)
nr <- as.numeric(system('wc -l cleaned_v3a.methylation | awk \'{print $1}\' ',intern=T))

#get list of numeric column classes (reads file faster)
dat_temp <- read.table("cleaned_v3a.methylation", header=T,nrows=50)
classes_to_use <- sapply(dat_temp , class)

#read file
dat0 <- read.table("cleaned_v3a.methylation", header=T,nrows=nr,colClasses=classes_to_use)


######## we want to filter out the cross hybridizing and SNP probes #######

############# load in the other stuff like cross hybridizing probes or whatever #########
snp_probes <- read.table('PCA/cpgs_within_10bp_of_SNP.table', header=F,stringsAsFactors=F)
names(snp_probes)[1] <- "CpG"

cross_hybridizing <- read.table('PCA/cpg-non-specific-probes-Illumina450k.txt',header=T,stringsAsFactors=F)
names(cross_hybridizing)[1] <- "CpG"

badprobes <- c(snp_probes$CpG, cross_hybridizing$CpG)

#also logit 2 transform the data to help make valid inferences 

dat <- logit2(dat0[-which(row.names(dat0) %in% badprobes),])



#make a dataframe of subject names that goes in the order they appear in the file
dat_order <- as.data.frame(names(dat))

#name the subjects the same as found in the matching column in the batch file
names(dat_order) <- c('methylome_id_fixed')

#establish a note of what the ordering is
dat_order$order <- c(1:(dim(dat)[2]))

batch <- read.csv("data_analysis/demo2_wholedata_methylation_x3c.csv", header=T,na.strings=c("NA","#N/A","N/A"))

#merge batch data with subject name order data
covariates_ordered0 <- merge(batch,dat_order,by='methylome_id_fixed')


#sort the data by order so the phenotype now lines up
covariates_ordered <- covariates_ordered0[order(covariates_ordered0$order),]

#read in studymax covariate data, to get the studymax ids (if doing studymax)
covars_studymax <- read.csv('data_analysis/demo2_studymaxnov0_methylation_x1.csv',header=T,na.strings=c("NA","#N/A","N/A"),stringsAsFactors=F)


##### load in probe annotations #####
probe_annotations <- read.csv('misc/humanmethylation450_15017482_v1-2.csv',header=T,skip=7,stringsAsFactors=F)

##### order this file by the probe order in the data #####
all_qced_probes <- data.frame(row.names(dat))
names(all_qced_probes)[1] <- "IlmnID"
all_qced_probes$order <- 1:dim(all_qced_probes)[1]
merged_probe_types0 <- merge(all_qced_probes,probe_annotations,by="IlmnID",all.x=TRUE)
probe_annotations_use <- merged_probe_types0[order(merged_probe_types0$order),]

n_shelf <- probe_annotations_use[probe_annotations_use$Relation_to_UCSC_CpG_Island == "N_Shelf",]$IlmnID
n_shore <- probe_annotations_use[probe_annotations_use$Relation_to_UCSC_CpG_Island == "N_Shore",]$IlmnID
island <- probe_annotations_use[probe_annotations_use$Relation_to_UCSC_CpG_Island == "Island",]$IlmnID
opensea <- probe_annotations_use[is.na(probe_annotations_use$Relation_to_UCSC_CpG_Island),]$IlmnID
s_shelf <- probe_annotations_use[probe_annotations_use$Relation_to_UCSC_CpG_Island == "S_Shelf",]$IlmnID
s_shore <- probe_annotations_use[probe_annotations_use$Relation_to_UCSC_CpG_Island == "S_Shore",]$IlmnID

#the gene grouping stuff depends on what gene it annotates to, i'll just take the first piece
unlist_split <- function(x, ...)
{
	toret <- unlist(strsplit(x, ...) )[1]
	return(t(toret))
}


probe_annotations_use$UCSC_RefGene_Group_simple <-  sapply(probe_annotations_use$UCSC_RefGene_Group, unlist_split, split = c(";"))

f_utr <- probe_annotations_use[probe_annotations_use$UCSC_RefGene_Group_simple == "5'UTR",]$IlmnID
t_utr <- probe_annotations_use[probe_annotations_use$UCSC_RefGene_Group_simple == "3'UTR",]$IlmnID
tss1500 <-  probe_annotations_use[probe_annotations_use$UCSC_RefGene_Group_simple == "TSS1500",]$IlmnID
tss200 <-  probe_annotations_use[probe_annotations_use$UCSC_RefGene_Group_simple == "TSS200",]$IlmnID
body  <-  probe_annotations_use[probe_annotations_use$UCSC_RefGene_Group_simple == "Body",]$IlmnID
exon  <-  probe_annotations_use[probe_annotations_use$UCSC_RefGene_Group_simple == "1stExon",]$IlmnID
intergenic  <-  probe_annotations_use[is.na(probe_annotations_use$UCSC_RefGene_Group_simple),]$IlmnID


############### get average methylation for each subject at each location ############
covariates_ordered$mean_overall <- apply(dat,2,mean,na.rm=T)

#gene features means
covariates_ordered$f_utr <- apply(dat[f_utr,],2,mean,na.rm=T)
covariates_ordered$t_utr <- apply(dat[t_utr,],2,mean,na.rm=T)
covariates_ordered$tss1500 <- apply(dat[tss1500,],2,mean,na.rm=T)
covariates_ordered$tss200 <- apply(dat[tss200,],2,mean,na.rm=T)
covariates_ordered$body <- apply(dat[body,],2,mean,na.rm=T)
covariates_ordered$exon <- apply(dat[exon,],2,mean,na.rm=T)
covariates_ordered$intergenic <- apply(dat[intergenic,],2,mean,na.rm=T)

#genome regions means
covariates_ordered$n_shelf <- apply(dat[n_shelf,],2,mean,na.rm=T)
covariates_ordered$n_shore <- apply(dat[n_shore,],2,mean,na.rm=T)
covariates_ordered$island  <- apply(dat[island,], 2,mean,na.rm=T)
covariates_ordered$opensea <- apply(dat[opensea,],2,mean,na.rm=T)
covariates_ordered$s_shelf <- apply(dat[s_shelf,],2,mean,na.rm=T)
covariates_ordered$s_shore <- apply(dat[s_shore,],2,mean,na.rm=T)




###### take differences between cases and controls in refgene groups and by cpg island relation

######## at v3
covariates_ordered_subset <- subset(covariates_ordered, visit == 2 & (PTSDbroad == 1 | PTSDbroad_highest == 0))
dat_subset <- dat[, names(dat) %in% covariates_ordered_subset$methylome_id_fixed]

summary(
	lm(mean_overall ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(f_utr ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(t_utr ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(tss1500 ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(tss200 ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(body ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(exon ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(intergenic ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)


summary(
	lm(n_shore ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(n_shelf ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(island ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(opensea ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(s_shelf ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)
summary(
	lm(s_shore ~ PTSDbroad + CD8T+CD4T+NK+Bcell+Mono + PC1_HGDP + PC2_HGDP + PC3_HGDP, data= covariates_ordered_subset)
)


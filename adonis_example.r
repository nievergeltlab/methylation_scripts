#Classify cases using adonis

#most of these are unneeded
library('CpGassoc')
library(minfi)
library(vegan)
library('CpGassoc')
library(superpc)
library(isva)
library(cluster)



#get number of rows in file (reads file faster)
nr <- as.numeric(system('wc -l cleaned_v3a.methylation | awk \'{print $1}\' ',intern=T))

#get list of numeric column classes (reads file faster)
dat_temp <- read.table("cleaned_v3a.methylation", header=T,nrows=50)
classes_to_use <- sapply(dat_temp , class)

#read file
dat <- read.table("cleaned_v3a.methylation", header=T,nrows=nr,colClasses=classes_to_use)


#make a dataframe of subject names that goes in the order they appear in the file
dat_order <- as.data.frame(names(dat))

#name the subjects the same as found in the matching column in the batch file
names(dat_order) <- c('methylome_id_fixed')

#establish a note of what the ordering is
dat_order$order <- c(1:(dim(dat)[2]))


batch <- read.csv("data_analysis/demo2_wholedata_methylation_x3i.csv", header=T,na.strings=c("NA","#N/A","N/A"))

#merge batch data with subject name order data
covariates_ordered0 <- merge(batch,dat_order,by='methylome_id_fixed')

#sort the data by order so the phenotype now lines up
covariates_ordered <- covariates_ordered0[order(covariates_ordered0$order),]

covars_ordered_use <- subset(covariates_ordered, visit == 3)


##### impute missing data to row median #####
impute_median <- function(x)
{
	med <- median(x,na.rm=T)
	x[is.na(x)] <- med
	return(x)
}

dat_subset <- dat[, names(dat) %in% covars_ordered_use$methylome_id_fixed]

dat_subset <- apply(dat_subset, 1, impute_median)
dat_subset <- logit2(dat_subset)




adamdonis <- adonis(dat_subset ~ DRRI_composite_i + CD4T+ CD8T + Bcell + Mono + NK + PC1_HGDP + PC2_HGDP + PC3_HGDP, data=covars_ordered_use, method="euclidean")


adamdonis3 <- adonis(dat_subset ~CD4T+ CD8T + Bcell + Mono + NK + PC1_HGDP + PC2_HGDP + PC3_HGDP +  DRRI_composite_i , data=covars_ordered_use, method="euclidean")


#I start off with an error 
#Error in vegdist(lhs, method = method, ...) :
#  NA/NaN/Inf in foreign function call (arg 1)

#vegdist is the distance calculation function. 
adonis(dat) ~ CD4T, data=batch)


#Note: data needs to be transposed to conform to NxP format.
#Note: seems to not take a matrix with missing data

adonis(dat_subset ~  CAPStots_v0 +  as.factor(PTSDbroad_highest) + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covars_ordered_use,method="euclidean" )


#From Schork:
#If one has collected N total samples and assayed the expression level of G genes on those samples, 
#then an N × N similarity matrix can be formed that reflects the correlation or similarity of the 
#samples with respect to the expression values over the G genes. 

#What kind of (Dis)similarity matrix do I make?
#Open question: How similar are any two nonrelated observations using any give metric?

#Comparison of daisy to vegdist on euclidean distance. 
#Conclusion: Program outputs match
#vegdist does not have a centering option
library(cluster)


testmat <- vegdist(dat_subset, method="euclidean")
testmatdaisy <- daisy(dat_subset)

#################################################################
##### Load all data


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

dat <- dat0[-which(row.names(dat0) %in% badprobes),]

#make a dataframe of subject names that goes in the order they appear in the file
dat_order <- as.data.frame(names(dat))

#name the subjects the same as found in the matching column in the batch file
names(dat_order) <- c('methylome_id_fixed')

#establish a note of what the ordering is
dat_order$order <- c(1:(dim(dat)[2]))

batch <- read.csv("data_analysis/demo2_wholedata_methylation_x3e.csv", header=T,na.strings=c("NA","#N/A","N/A"))

#merge batch data with subject name order data
covariates_ordered0 <- merge(batch,dat_order,by='methylome_id_fixed')

#sort the data by order so the phenotype now lines up
covariates_ordered <- covariates_ordered0[order(covariates_ordered0$order),]


#################################################################
##### Comparison of post deployment studymax data

#Filter data and impute missing values

covariates_ordered_subset <- subset(covariates_ordered, visit == 3 & (PTSDbroad == 1 | PTSDbroad_highest == 0))
dat_subset0 <- dat[, names(dat) %in% covariates_ordered_subset$methylome_id_fixed]

##### impute missing data to row median #####
impute_median <- function(x)
{
	med <- median(x,na.rm=T)
	x[is.na(x)] <- med
	return(x)
}

dat_subset <- apply(dat_subset0, 1, impute_median)
dat_subset <- logit2(dat_subset)

#Calculate euclidean distance
 dat_dist  <- vegdist(dat_subset, method="euclidean")

#No covariates
 adonis(dat_dist ~ PTSDbroad, data=covariates_ordered_subset, method="euclidean")
#All covariates
 adonis(dat_dist ~  PTSDbroad + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")

#Change covariate order each time
 adonis(dat_dist ~  PTSDbroad + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")
 adonis(dat_dist ~  DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP+ PTSDbroad  , data=covariates_ordered_subset, method="euclidean")
 adonis(dat_dist ~   age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP+ PTSDbroad + DRRI_composite_i   , data=covariates_ordered_subset, method="euclidean")
 adonis(dat_dist ~  CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP+ PTSDbroad + DRRI_composite_i +  age_v0   , data=covariates_ordered_subset, method="euclidean")
 adonis(dat_dist ~  CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP+ PTSDbroad + DRRI_composite_i +  age_v0 + CD8T  , data=covariates_ordered_subset, method="euclidean")
 adonis(dat_dist ~  NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP+ PTSDbroad + DRRI_composite_i +  age_v0 + CD8T + CD4T , data=covariates_ordered_subset, method="euclidean")
 adonis(dat_dist ~  Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP+ PTSDbroad + DRRI_composite_i +  age_v0 + CD8T + CD4T + NK, data=covariates_ordered_subset, method="euclidean")
 adonis(dat_dist ~  Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP+ PTSDbroad + DRRI_composite_i +  age_v0 + CD8T + CD4T + NK + Bcell, data=covariates_ordered_subset, method="euclidean")
 adonis(dat_dist ~  PC1_HGDP+PC2_HGDP+PC3_HGDP+ PTSDbroad + DRRI_composite_i +  age_v0 + CD8T + CD4T + NK + Bcell + Mono , data=covariates_ordered_subset, method="euclidean")
 adonis(dat_dist ~  PC2_HGDP+PC3_HGDP+ PTSDbroad + DRRI_composite_i +  age_v0 + CD8T + CD4T + NK + Bcell + Mono + PC1_HGDP , data=covariates_ordered_subset, method="euclidean")
 adonis(dat_dist ~  PC3_HGDP + PTSDbroad + DRRI_composite_i +  age_v0 + CD8T + CD4T + NK + Bcell + Mono + PC1_HGDP + PC2_HGDP, data=covariates_ordered_subset, method="euclidean")

#CTQ
 adonis(dat_dist ~  CTQ_TOTAL + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")


 #reordered
  posttest <- adonis(dat_dist ~  CAPStots_v0 +   DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP+as.factor(PTSDbroad_highest) , data=covariates_ordered_subset, method="euclidean")

#No cell types
 posttest2 <- adonis(dat_dist ~  CAPStots_v0 +  as.factor(PTSDbroad_highest) + DRRI_composite_i + age_v0+ PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")
#No v0 pTSD
 posttest2 <- adonis(dat_dist ~  as.factor(PTSDbroad_highest) + DRRI_composite_i + age_v0+ PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")
#Other PTSD phenotype
 posttest <- adonis(dat_dist ~  CAPStots_v0 +  as.factor(CAPSdx_highest) + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")
#tobac
 posttest <- adonis(dat_dist ~  CAPStots_v0 +  TOBAC_amtuse_est_i + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")
 #tobac minus covars
  adonis(dat_dist ~  TOBAC_amtuse_est_i + DRRI_composite_i + age_v0+PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")

#alc
 adonis(dat_dist ~  CAPStots_v0 + AUDIT_combined_i + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")




#take only the most variable positions
 probesds <- apply(dat_subset,2,sd)
 probesds2 <- sort(probesds,decreasing=TRUE)
 deviant_probes <- names(probesds2)[1:(dim(dat_subset)[2]*.66)]

 dat_dist_deviant <- vegdist(dat_subset[,colnames(dat_subset) %in% deviant_probes], method="euclidean")

 deviant_probes2 <- names(probesds2)[1:(dim(dat_subset)[2]*.10)]
 dat_dist_deviant2 <- vegdist(dat_subset[,colnames(dat_subset) %in% deviant_probes], method="euclidean")


 posttest_deviant <- adonis(dat_dist_deviant ~   as.factor(PTSDbroad_highest) + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")
 posttest_deviant2 <- adonis(dat_dist_deviant2 ~  CAPStots_v0 +  as.factor(PTSDbroad_highest) + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")

#centered euclidean distance data
 dat_dist_center  <- daisy(dat_subset, metric="euclidean", stand=TRUE)
  adonis(dat_dist_center ~  PTSDbroad_highest + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")


#jaccard distance data
 #results may be meaningless because data have negative entries in method âjaccardâ

 dat_dist_jaccard  <- vegdist(dat_subset, method="ja")


 posttest_jaccard <- adonis(dat_dist_jaccard ~  CAPStots_v0 +  as.factor(PTSDbroad_highest) + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")


#kulcuzinski
 #results may be meaningless because data have negative entries in method âjaccardâ
 dat_dist_ku <- vegdist(dat_subset, method="ku")

 posttest_ku <- adonis(dat_dist_ku ~  CAPStots_v0 +  as.factor(PTSDbroad_highest) + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="ku")

#Different metric on B values
 dat_dist_b  <- vegdist(ilogit2(dat_subset), method="euclidean")

#All covariates
 adonis(dat_dist_b ~  CAPStots_v0 +  as.factor(PTSDbroad_highest) + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="euclidean")

 dat_dist_ku <- vegdist(ilogit2(dat_subset), method="ku")

 adonis(dat_dist_ku~  CAPStots_v0 +  as.factor(PTSDbroad_highest) + DRRI_composite_i + age_v0+ CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=covariates_ordered_subset, method="ku")


#Filter data and impute missing values

covariates_ordered_subset <- subset(covariates_ordered, visit == 0)
dat_subset0 <- dat[, names(dat) %in% covariates_ordered_subset$methylome_id_fixed]

##### impute missing data to row median #####
impute_median <- function(x)
{
	med <- median(x,na.rm=T)
	x[is.na(x)] <- med
	return(x)
}

dat_subset <- apply(dat_subset0, 1, impute_median)
dat_subset <- logit2(dat_subset)



#Calculate euclidean distance
 dat_dist  <- vegdist(dat_subset, method="euclidean")

#All covariates
 adonis(dat_dist ~   CD8T+CD4T+NK+Bcell+Mono +PC1_HGDP+PC2_HGDP+PC3_HGDP + PTSDbroad_highest +  age_v0, data=covariates_ordered_subset, method="euclidean")


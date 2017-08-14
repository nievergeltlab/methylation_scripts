#Classify methylation status using unsupervised PCA

library(glmnet)
library('CpGassoc')
library(minfi)
library(superpc)
library(isva)

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

batch0 <- read.csv("data_analysis/demo2_wholedata_methylation_x3d.csv", header=T,na.strings=c("NA","#N/A","N/A"))


ages <- read.csv('DNAm/methylation_swan_noqc_nodnambmiq_6_29_2015.csv',header=T)
names(ages)[1] <- "methylome_id_fixed"

batch <- merge(ages, batch,by="methylome_id_fixed")
batch$agedif <- batch$age_i - batch$DNAmAge



#merge batch data with subject name order data
covariates_ordered0 <- merge(batch,dat_order,by='methylome_id_fixed')

#sort the data by order so the phenotype now lines up
covariates_ordered <- covariates_ordered0[order(covariates_ordered0$order),]

#read in studymax covariate data, to get the studymax ids (if doing studymax)
covars_studymax <- read.csv('data_analysis/demo2_studymaxnov0_methylation_x1.csv',header=T,na.strings=c("NA","#N/A","N/A"),stringsAsFactors=F)

######################################################################
######################################################################
# 6 mo analysis in all subjects
######################################################################
######################################################################
overlapping_27k <- read.csv('DNAm/datMiniAnnotation27k.csv', header=T)

covariates_ordered_subset <- subset(covariates_ordered) #, visit == 3 & (PTSDbroad == 1 | PTSDbroad_highest == 0))
dat_subset0 <- dat[row.names(dat) %in% overlapping_27k$Name, names(dat) %in% covariates_ordered_subset$methylome_id_fixed]
covs <- subset(covariates_ordered_subset, select=c(age_v0, CD8T,CD4T,NK,Bcell,Mono ,PC1_HGDP,PC2_HGDP,PC3_HGDP))

##### impute missing data to row median #####
impute_median <- function(x)
{
	med <- median(x,na.rm=T)
	x[is.na(x)] <- med
	return(x)
}

dat_subset <- t(apply(dat_subset0, 1, impute_median))
dat_subset <- logit2(dat_subset)



######################################################################
# Before residualizing data
######################################################################

#Get number of significant PCs

#rescale the data first, has to be transposed because this process is done via column
dat_rescale <- scale(t(dat_subset),center=TRUE,scale=TRUE)


sigpcs <- EstDimRMT(t(dat_rescale))
sigpcs$dim

# Do pcs on the data. Q mode pca (more variables than subjects)
result <- prcomp(t(dat_subset),scale=TRUE)

# variance explained
summary(result)

#extract pcs
PCs <- as.data.frame(result$x)
PCs1 <- PCs
PCs1$methylome_id_fixed <- row.names(PCs1)


pc_use <- merge(PCs1,covariates_ordered_subset,by="methylome_id_fixed")

#plot the PC
pc_use$color <- NA
pc_use$color <- abs(pc_use$agedif)/max(abs(pc_use$agedif))

outlier_ids <- scan(what="character")
4282_0
4803_3
4786_2
4275_3
4122_2
4393_2
4255_0
4817_3
4122_0
4347_3
4275_2
4311_3
4393_3

##Plot PCs
plot(pc_use$PC1, pc_use$PC2, col='white' )
print(
	text(pc_use$PC1,pc_use$PC2,labels=pc_use$id_visit,col=rgb(pc_use$color,0,pc_use$color,1),cex=1.5)
)
print(
	text(pc_use$PC1,pc_use$PC2,labels=pc_use$id_visit)
)

print(
	text(pc_use[pc_use$id_visit %in% outlier_ids,]$PC1,pc_use[pc_use$id_visit %in% outlier_ids,]$PC2,labels=pc_use[pc_use$id_visit %in% outlier_ids,]$id_visit,col="red")
)


#Test a PC for association with PTSD

summary(
	lm(PC18 ~ NK + Bcell + Mono + CD4T + CD8T + PTSDbroad + DRRI_composite_i+ age_v0 +PC1_HGDP+PC2_HGDP+PC3_HGDP, data=subset(pc_use, PC18 < 250))
)

#correlate each PC with PTSD
correlations_r <- corr.test(cbind( covariates_ordered_subset$PTSDbroad,PCs1[,1:26]),method="spearman")$r
correlations_p <- corr.test(cbind( covariates_ordered_subset$PTSDbroad,PCs1[,1:26]),method="spearman")$p

#Adjust for multiple comparisons
p.adjust(sort(correlations_p[-1,1]))

#Plot significant PCs
plot(pc_use$PC18, pc_use$PC14, col=pc_use$color )

#whichm methylation sites predict this PC?
pcpred <- cpg.assoc(dat, PCs1$PC18,  logit.transform = FALSE, chip.id = NULL, subset = NULL, random = FALSE, fdr.cuto = 0.05, large.data = TRUE, fdr.method = "BH", logitperm= FALSE)


######################################
######### correlate each PCA to probes
################################################

summary(
	glm(PC7 ~ NK, data=subset(pc_use),family=binomial)
)
dat_rescale2 <- data.frame(dat_rescale)

dat_rescale2$methylome_id_fixed <- row.names(dat_rescale2)

d_a <- merge(pc_use, dat_rescale2, by="methylome_id_fixed")

summary(
	lm(PC1 ~ PC1_HGDP +  cg03118626, data=d_a)
)
summary(
	lm(PC1 ~ as.factor(plate), data=d_a)
)

plot(d_a$PC7 , d_a$cg03118626)

#load the top 100 probes and remove them and recalculate PCA...
file=paste('data_analysis/RESULTS/all_6mo_condv0_PTSDbroad_noFC_agecellpcdrri.mwas')

#get list of numeric column classes (reads files faster)
dat_temp <- read.table(file, header=T,nrows=50,stringsAsFactors=F)
classes_to_use <- sapply(dat_temp , class)

#get number of rows for largest file
nr <- as.numeric(system('wc -l data_analysis/RESULTS/all_6mo_condv0_PTSDbroad_noFC_agecellpcdrri.mwas | awk \'{print $1}\' ',intern=T))

#read file and transpose
ewas_results <- read.table(file, header=F,colClasses=classes_to_use,nrows=nr)
names(ewas_results ) <- c("case","IlmnID","beta","se","z","p")
ewas_results <- ewas_results[order(ewas_results0$p),]
ewas_results[1:100,]$IlmnID




 



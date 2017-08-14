#Semi supervised RPMM models

source('SSRPMM.R')

#####################################################################################################
################################## Read in logit transformed methylation values #######################################
#####################################################################################################
nr <- as.numeric(system('wc -l cleaned_postcombat_v1b_cell_removed_logit.methylation | awk \'{print $1}\' ',intern=T))

#get list of numeric column classes (reads file faster)
dat_temp <- read.table("cleaned_postcombat_v1b_cell_removed_logit.methylation", header=T,nrows=50)
classes_to_use <- sapply(dat_temp , class)

#read file
dat <- read.table("cleaned_postcombat_v1b_cell_removed_logit.methylation", header=T,nrows=nr,colClasses=classes_to_use)

#transpose the data, and logit transform it (if not already)
dat_t <- t(dat)

#get a subject pool
covars <- read.csv('data_analysis/demo2_wholedata_methylation_x2.csv',header=T,na.strings=c("NA","#N/A","N/A"),stringsAsFactors=F)

#take the 6 months analysis because marco has that, do not include potential false controls
covars_use <- subset(covars, visit == 3 & (PTSDbroad == 1 | PTSDbroad_highest == 0))

dat_t_use <- dat_t[row.names(dat_t) %in% covars_use$methylome_id_fixed,]


#covariates may have to be ordered the same as methylation data so we do that now

sample_ordering <- data.frame(cbind(row.names(dat_t_use) , 1:length(row.names(dat_t_use))), stringsAsFactors=F)
names(sample_ordering) <- c("methylome_id_fixed", "ordering")
#for some reason the ordering is a string
sample_ordering$ordering <- as.numeric(sample_ordering$ordering)

#merge the sample ordering into the covariat edata 
case_order <- merge(sample_ordering, covars_use,by="methylome_id_fixed")
#sort the covariates by the sample ordering
case_order2 <- case_order[order(case_order$ordering),]

#make case status subset
case_status <- subset(case_order2, select="PTSDbroad")
#make a covariate subset
covariates_use <- subset(case_order2, select=c(PTSDbroad,age_v0,PC1_HGDP,PC2_HGDP,PC3_HGDP))



J1 = dim(dat_t_use)[2] #number of CpG loci in the HNSCC dataset
N1 = dim(dat_t_use)[1] #number of samples in the HNSCC dataset
P1 = dim(covariates_use)[2] #number of covariate factors in the HNSCC covariate data
N1 == dim(covariates_use)[1] #should be true!

#split the data randomly between training and test sets, using PTSD broad to do a stratified split.
split_dat <-  TrainTestSplit(dat_t_use, covariates_use, Strat= "PTSDbroad", seed= 17)

########extract training and testing sets
data_training <- split_dat[[1]]
data_testing <- split_dat[[2]]

#Training Data
data_trainingBetas = data_training[,-(1:(P1+1))]
data_trainingCovariates = data_training[,(1:(P1 + 1))]
#Testing Data
data_testingBetas = data_testing[,-(1:(P1+1))]
data_testingCovariates = data_testing[,(1:(P1 + 1))]

########## get important cpgs
data_trainingScores = MostImpCpGs(Y = data_trainingBetas, covariates = data_trainingCovariates,clinvar = "PTSDbroad", 
	terms = c("age_v0" ,"PC1_HGDP" , "PC1_HGDP" , "PC3_HGDP"), factors = NULL,is.factor=TRUE)

data_training_CXvalidationResults = NestedXValidation(Y = data_trainingBetas, covariates = data_trainingCovariates, 
	TScores = data_trainingScores, clinvar = "PTSDbroad", vartype = "binary", mrange = c(5,50), method = "gaussian", L = 20, seeds = 1:20)

#get M with lowest pv 
sort(apply(data_training_CXvalidationResults, 1, median,na.rm=T))

mrange = 5:20
loessCurve = loess.smooth(mrange, data_training_CXvalidationResults, degree = 2)
MOpt_data_training = subset(data.frame(loessCurve$x, loessCurve$y), loessCurve$y ==
min(loessCurve$y))[[1]]
par(mar = c(5,5,4,2))
plot(mrange, data_training_CXvalidationResults, cex = 0.75, xlab = "Number of top ranking loci (M)",
ylab = "Median P-value", cex.lab = 2, cex.axis = 1.5)
lines(loessCurve$x, loessCurve$y, lwd = 5)
abline(v = MOpt_data_training, col = "red", lwd = 2, lty = "dashed")


data_ClassesTesting = PredMethClasses(Ytrain = data_trainingBetas,Ytest = data_testingBetas, Scores = data_trainingScores, M = 10,
method = "gaussian")

permTestChiSquare(data_ClassesTesting, data_testingCovariates[,"PTSDbroad"])


#Determine differences in sample profiles using Euclidean Distance

library(plyr)

#for each subject, get a distance pairing
#residualize each first?

#Should data be rescaled together or per visit?

#take only subjects with V0 and V3 data
usable <- intersect(subset(covariates_ordered, visit == 0)$studyid,subset(covariates_ordered, visit == 3)$studyid)

covariates_ordered_subset <- subset(covariates_ordered, (visit == 3 | visit == 0) & studyid %in% usable )
dat_subset0 <- dat[, names(dat) %in% covariates_ordered_subset$methylome_id_fixed]

##### impute missing data to row median #####
impute_median <- function(x)
{
	med <- median(x,na.rm=T)
	x[is.na(x)] <- med
	return(x)
}

dat_subset <- t(apply(dat_subset0, 1, impute_median))
dat_subset <- logit2(dat_subset)

dat_subset2 <- cbind(as.data.frame(covariates_ordered_subset$studyid),t(dat_subset))
names(dat_subset2)[1] <- "studyid"



distfun  <- function(covar_df, df)
{
	use <- df[names(df) %in% covar_df$methylome_id_fixed,]
	return(dist(rbind(use[,1],use[,2])))
}

euclidean_distances <- ddply(covariates_ordered_subset, ~ studyid, distfun,df=dat)
names(euclidean_distances)[1] <- "studyid"

names(euclidean_distances)[2] <- "distance"

covs_v3 <- subset(covariates_ordered_subset, visit == 3)

d1 <- merge(euclidean_distances,covs_v3,by="studyid")

summary(
	lm(distance ~ PTSDbroad,data=d1)
)

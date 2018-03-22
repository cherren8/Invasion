# Analysis to evaluate whether relapsing FMT patients differ from cured patients in gut community connectivity 
# CMH 23Jan18

####################create necessary functions######################

#find the number of zeroes in a vector
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}

#create function that averages only negative values in a vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}

#create function that averages only positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}

min.non.zero <- function(vec) {
  smallest <- min(vec[vec != 0])
  return(smallest)
}

###################################################################
###################################################################
### Workflow options ####
###################################################################
###################################################################

## Choose a persistence cutoff (min. fraction of taxon presence) for retaining taxa in the analysis
pers.cutoff <- 0.6 #  .2, .3, .4, .6 with .001 showed significant pos/neg. use .4, though .6 gives best classification 
## Choose a mean value cutoff
mean.cutoff <- 0.0001 #0.001 
## Decide the number of iterations to run for each taxon. (>= 200 is recommended)
# Larger values of iter mean the script takes longer to run
iter <- 200
## Decide whether to use taxon/column shuffle (col.shuffle = T) or row shuffle algorithm (col.shuffle = F)
col.shuffle <- F
#Calculate cohesion scores for relapsing and cured patients
checkup.day <- 28
#decide whether to evaluate donor cohesion or patient cohesion
test.donor.coh <- F
#Decide whether to include all donors or only donors of analyzed patients
use.only.an <- T #T or F similar. Have been using T 

###################################################################
###################################################################

# Read in dataset
## Data should be in a matrix where each row is a sample. 
b.orig <- read.csv("recurring_cdiff_09Oct17.csv", header = T, row.names = 1) #138, 4643 closed ref OTU picking, no deblurring 

dim(b.orig)

#set universal detection limit across all samples to account for different sequencing depth 
b.orig <- b.orig[rowSums(b.orig) > 0, colSums(b.orig) > 0 ]
min.sens.all <- max(apply(b.orig , 2, min.non.zero))
b.orig[b.orig < min.sens.all] <- 0
dim(b.orig)

meta <- read.csv("RecurringCdiffMeta_08Oct17.csv", header = T, row.names = 1)


#Subset to only donors
donor.samp <- paste0("X", rownames(meta)[meta$donor_or_patient == "Donor"])
colnames(b.orig)
b <- b.orig[, colnames(b.orig) %in% donor.samp]
dim(b)

#optionally, subset to only donors of analyzed patients
if(use.only.an){
  exp.code.an <- unique(meta$exp_code[meta$day_since_fmt == 28]) #pulls experiment code for patients at day 28
  b <- b[, meta$exp_code[match(colnames(b), paste0("X", rownames(meta)) )] %in% exp.code.an ]
  dim(b)
}

# Suggested steps to re-format data. At the end of these steps, the data should be in a matrix "c" where there are no empty samples or blank taxon columns. 
c <- t(as.matrix(b))
c <- c[rowSums(c) > 0, colSums(c) > 0]
dim(c)

# Save total number of individuals in each sample in the original matrix.
rowsums.orig <- rowSums(c)

# Based on persistence cutoff, define a cutoff for the number of zeroes allowed in a taxon's distribution
zero.cutoff <- ceiling(pers.cutoff * dim(c)[1])
  
# Remove taxa that are below the persistence cutoff
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
# Remove any samples that no longer have any individuals, due to removing taxa
d <- d[rowSums(d) > 0, ]

# Create relative abundance matrix.  
rel.d <- d / rowsums.orig

#Remove low abundance taxa
rel.d <- rel.d[, colMeans(rel.d) > mean.cutoff]

# Optionally, check to see what proportion of the community is retained after cutting out taxa
par(mfrow = c(1, 1))
hist(rowSums(rel.d))

# Create observed correlation matrix
cor.mat.true <- cor(rel.d)

# Create vector to hold median otu-otu correlations for initial otu
med.tax.cors <- vector()
dim(rel.d)

set.seed(88)
seeds <- sample(seq(1:100000000), size = (dim(rel.d)[2] * iter * dim(rel.d)[1]))

# Run this loop for the null model to get expected pairwise correlations
if(col.shuffle) {
  for(which.taxon in 1:dim(rel.d)[2]){
    
    #create vector to hold correlations from every permutation for each single otu
    ## perm.cor.vec.mat stands for permuted correlations vector matrix
    perm.cor.vec.mat <- vector()
    
    for(i in 1:iter){
      #Create empty matrix of same dimension as rel.d
      perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
      rownames(perm.rel.d) <- rownames(rel.d)
      colnames(perm.rel.d) <- colnames(rel.d)
      
      #For each otu
      for(j in 1:dim(rel.d)[2]){ 
        #Set new seed for upcoming randomization
        set.seed(seeds[j + (i-1)*(dim(rel.d)[1]) + (which.taxon-1)*dim(rel.d)[2] ])
        # Replace the original taxon vector with a permuted taxon vector
        perm.rel.d[, j ] <- sample(rel.d[ ,j ]) 
      }
      
      # Do not randomize focal column 
      perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]
      
      # Calculate correlation matrix of permuted matrix
      cor.mat.null <- cor(perm.rel.d)
      
      # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
      perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
      
    }
    # Save the median correlations between the focal taxon and all other taxa  
    med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
    
    # For large datasets, this can be helpful to know how long this loop will run
    if(which.taxon %% 20 == 0){print(which.taxon)}
  }
} else {
  for(which.taxon in 1:dim(rel.d)[2]){
    
    #create vector to hold correlations from every permutation for each single otu
    ## perm.cor.vec.mat stands for permuted correlations vector matrix
    perm.cor.vec.mat <- vector()
    
    for(i in 1:iter){
      #Create duplicate matrix to shuffle abundances
      perm.rel.d <- rel.d 
      
      #For each taxon
      for(j in 1:dim(rel.d)[1]){ 
        which.replace <- which(rel.d[j, ] > 0 ) 
        # if the focal taxon is greater than zero, take it out of the replacement vector, so the focal abundance stays the same
        which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
        #Set new seed for upcoming randomization
        set.seed(seeds[j + (i-1)*(dim(rel.d)[1]) + (which.taxon-1)*dim(rel.d)[2] ])
        #Replace the original taxon vector with a vector where the values greater than 0 have been randomly permuted 
        perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal]) 
      }

      # Calculate correlation matrix of permuted matrix
      cor.mat.null <- cor(perm.rel.d)
      
      # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
      perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
      
    }
    # Save the median correlations between the focal taxon and all other taxa  
    med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
    
    # For large datasets, this can be helpful to know how long this loop will run
    if(which.taxon %% 20 == 0){print(which.taxon)}
  }
 }

  
# Save observed minus expected correlations. Use custom correlations if use.custom.cors = TRUE
obs.exp.cors.mat <- cor.mat.true - med.tax.cors
diag(obs.exp.cors.mat) <- 0

#### 
#### Produce desired vectors of connectedness 

# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

par(mfrow = c(1, 2))
hist(connectedness.neg)
hist(connectedness.pos)

#Create connectedness vectors for all OTUs by inserting zeroes 
conn.pos.all <- rep(0, dim(b.orig)[1])
conn.pos.all[match(names(connectedness.pos), rownames(b.orig))] <- connectedness.pos

conn.neg.all <- rep(0, dim(b.orig)[1])
conn.neg.all[match(names(connectedness.neg), rownames(b.orig))] <- connectedness.neg

#find patient IDs to calculate whether donor cohesion influences success
relapsePID <- paste0("X", rownames(meta)[meta$number_recurrence_after_fmt == 1 & meta$day_since_fmt == checkup.day])
curedPID <- paste0("X", rownames(meta)[meta$number_recurrence_after_fmt != 1 & meta$day_since_fmt == checkup.day])

#find donor IDs to calculate whether donor cohesion influences success
relapse.donors <- meta$exp_code[rownames(meta) %in% rownames(meta)[meta$number_recurrence_after_fmt == 1 & meta$day_since_fmt == checkup.day]] 
cure.donors <- meta$exp_code[rownames(meta) %in% rownames(meta)[meta$number_recurrence_after_fmt != 1 & meta$day_since_fmt == checkup.day]] 
relapseDID <- paste0("X", rownames(meta)[meta$donor_or_patient == "Donor" & meta$exp_code %in% relapse.donors])
curedDID <- paste0("X", rownames(meta)[meta$donor_or_patient == "Donor" & meta$exp_code %in% cure.donors])

#Calculate cohesion at day 28 for cured and relapsing

ifelse(test.donor.coh, {
  curedID <- curedDID
  relapseID <- relapseDID
}, {curedID <- curedPID
    relapseID <- relapsePID} 
)
  
b.cure <- t(b.orig[, colnames(b.orig) %in% curedID])
dim(b.cure)
b.relapse <- t(b.orig[, colnames(b.orig) %in% relapseID])
dim(b.relapse)

#Combine cured and relapsing patients into one matrix
b.both <- rbind(b.cure, b.relapse)

coh.neg.both <- b.both %*% conn.neg.all
coh.pos.both <- b.both %*% conn.pos.all

cr.factor <- c(rep("cure", times = dim(b.cure)[1]), rep("relapse", times= dim(b.relapse)[1]))

par(mfrow = c(1, 1))
plot(coh.pos.both, coh.neg.both, col = as.factor(cr.factor))
 
#Calculate whether difference from host predicts relapse
#Need to record as NA if donor sample not present 
id.both <- c(curedID, relapseID)
pair.code <- meta$donor_group[match(id.both , paste0("X", rownames(meta)) )]
length(pair.code) #Corresponds to which donor the patient received
exp.code <- meta$exp_code[match(id.both , paste0("X", rownames(meta)) )] #unique code for each patient 

bc.pre.28 <- vector()

for(z in 1:length(pair.code)){
  i <- pair.code[z]
  
  donor.temp.id <- rownames(meta)[meta$donor_group == i & meta$donor_or_patient == "Donor" & meta$exp_code == exp.code[z]]
  patient.temp.id <- rownames(meta)[meta$donor_group == i & meta$donor_or_patient == "Patient" & meta$day_since_fmt == checkup.day & meta$exp_code == exp.code[z]]
  
  donor.samp.temp <- b.orig[, colnames(b.orig) == paste0("X", donor.temp.id)]
  patient.samp.temp <- b.orig[, colnames(b.orig) == paste0("X", patient.temp.id)]
  
  if(dim(as.matrix(donor.samp.temp))[2] > 0 & dim(as.matrix(patient.samp.temp))[2] > 0) { 
    bc.pre.28[z] <- vegdist(rbind(donor.samp.temp, patient.samp.temp))
  } else {
      bc.pre.28[z] <- NA }
}

summary(lm(bc.pre.28 ~ cr.factor))


par(mfrow = c(1, 2))
par(bty = "l")
boxplot(coh.pos.both ~ cr.factor, ylab = "Positive Cohesion")
stripchart(coh.pos.both ~ cr.factor, add = T, vertical = T)
boxplot(coh.neg.both ~ cr.factor, ylab = "Negative Cohesion", bty = "l")
stripchart(coh.neg.both ~ cr.factor, add = T, vertical = T)


b.sub <- t(b.orig[, match(id.both, colnames(b.orig))])
dim(b.sub)
b.sub <- b.sub[, match(names(connectedness.pos), colnames(b.sub))]
dim(b.sub)
dim(rel.d)


library(vegan)
#Check to make sure diversity is not driven by depth variation. Calculate diversity after imposing universal detection limit
b.rare <- b.both
min.sens <- max(apply(b.rare, 1, min.non.zero))
b.rare[b.rare < min.sens] <- 0
div.both <- apply(b.rare, 1, diversity)

summary(lm(bc.pre.28 ~ cr.factor))
summary(lm(div.both ~ cr.factor))

summary(lm(coh.pos.both ~ cr.factor)) #Significant for parametric ANOVA, also
summary(lm(coh.neg.both ~ cr.factor))

wilcox.test(coh.pos.both ~ cr.factor , alternative = "greater", exact = T)
wilcox.test(coh.neg.both ~ cr.factor , alternative = "less", exact = T)



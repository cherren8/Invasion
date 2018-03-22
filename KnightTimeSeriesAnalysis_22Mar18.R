## Knight Lab Daily Stool Time Series
## 30 Jan 18 

library(lme4)
library(merTools)
library(MuMIn)
library(lmerTest)

########create necessary functions ##########
select.vec <- function(n, vec){
  back <- vec[n]
  return(back)
}

substring.last <- function(string){
  length <- nchar(string)
  last <- substring(string, first = length, last = length)
  return(last)
}

zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}

#create function that averages only negative values in a vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- sum(neg.vals) / length(vector)
  return(n.mean)
}

#create function that averages only positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0 & vector < 1)]
  p.mean <- sum(pos.vals) / length(vector)
  return(p.mean)
}

prop.over <- function(vec, cutoff){
  amt.over <- length(which(vec > cutoff)) / length(vec)
  return(amt.over)
}

#Create a function to calculate longest streak in each OTU
streakval <- function(vec){
  length.streak <- rep(0, length(vec))
  
  for(i in 1:length(vec)){
    j <- i
    while(vec[j] > 0 & j <= length(vec)){
      length.streak[i] <- length.streak[i] + 1
      j <- j + 1
    }
  }
  return(max(length.streak))
}

#Create a function to calculate probability of longest streak in each OTU
streakp <- function(vec){
  p.var <- mean(ceiling(vec) > 0)
  q.var <- 1 - p.var
  r.var <- max(1, streakval(vec) )
  n.var <- length(vec)
  
  poly.roots <- signif(polyroot(c(1, -1, rep(0, r.var - 1), q.var*(p.var^r.var)) ), digits = 15)
  poly.roots.real <- Re(poly.roots)[abs(Im(poly.roots)) < 10^-10]
  
  x.var <- poly.roots.real[abs(poly.roots.real - 1) == min(abs(poly.roots.real - 1))]
  
  streak.pval <- 1 - min( ((1 - p.var * x.var) /( (r.var + 1 - r.var*x.var )*q.var ) ) * (1 / (x.var^ (n.var + 1))), 1)
  
  return(streak.pval)
}

############################################################################################# 
############################################################################################# 

#must use the version of the data that matches the taxonomy file
trim <- read.csv("ClosedRef_10Dec17.csv", header = T, row.names = 1)

meta <- read.csv("AllSeriesMeta_09Dec17.csv", header = T, row.names = 1)

dim(trim)
#Change OTU names to integers
ab <- trim
colnames(ab) <- paste(rep("otu", times = dim(ab)[2]), seq(1, dim(ab)[2], 1), sep = ".")
dim(ab)
colnames(ab)[1:5]
rownames(ab)[1:5]
rownames(meta)[1:5]

#Put samples in same order as metadata
meta.in <- meta[rownames(meta) %in% rownames(ab), ]
dim(meta.in)

as.Date(meta.in$collection_timestamp[1:20], format = "%m/%d/%Y")
meta.ord <- meta.in[order(as.Date(meta.in$collection_timestamp, format = "%m/%d/%Y")), ]

ab.ord <- ab[match(rownames(meta.ord), rownames(ab)), ]
mean(rownames(ab.ord) == rownames(meta.ord)) #rownames are equal

#Subset to stool samples from subject M03 
m <- ab.ord[meta.ord$host_subject_id == "M03" & meta.ord$env_material == "feces", ]
m.meta <- meta.ord[meta.ord$host_subject_id == "M03" & meta.ord$env_material == "feces", ]
m.meta$date <- as.Date(m.meta$collection_timestamp, format = "%m/%d/%Y")

#Subset to years 2013 and later
m <- m[as.numeric(format(m.meta$date,'%Y')) >= 2013, ]
m.meta <- m.meta[as.numeric(format(m.meta$date,'%Y')) >= 2013, ]
dim(m)

plot(julian(m.meta$date))
#mostly continuous except for a couple gaps. Use the longest string of consecutive time points for analysis

m2 <- m[608:dim(m)[1], ]
m2.meta <- m.meta[608:dim(m)[1], ]
summary(m2.meta$geo_loc_name)
dim(m2)

############################### input options for analysis
#pick a minimum proportion of samples to analyze OTU
#Significant results for m2 are .4, F, .0001, 1, F, -9
zero.prop <- .5
use.prop.cutoff <- F
#choose mean abundance cutoff
mean.cutoff <- 0.0001 #.0001 has better fit via AIC than .001 for both .01 and .05 streak cutoff
#define the time lag 
time.lag <- 1 
col.shuffle <- F #(T also gives similar results)
det.cutoff <- -9

#########################################################

a <- m2

colnames.save <- colnames(a)
rownames.save <- rownames(a)
                
a.num <- as.matrix(a)
hist(rowSums(a)) #slight integer errors mean that they do not sum exactly to 1

a.num <- a.num[rowSums(a.num) > 0, colSums(a.num) > 0]

#Now, create global detection limit
det.mat <- a.num
det.mat[det.mat == 0] <- 1
hist(log(apply(det.mat, 1, min)), breaks = 50)

#Take out samples with detection lower than universal detection cutoff
a.mat <- a.num[log(apply(det.mat, 1, min)) < det.cutoff, ]
a.mat[a.mat < exp(det.cutoff)] <- 0
hist(rowSums(a.mat))
dim(a.mat)

b <- a.mat

#create new dataframe to manpiulate
d.sub <- as.matrix(a.mat)

#take out OTUs with fewer occurrences than the presence cutoff 
zero.cutoff <- ceiling(zero.prop*dim(d.sub)[1])
if(use.prop.cutoff) d.sub <- d.sub[, apply(d.sub, 2, zero) < (dim(d.sub)[1]-zero.cutoff) ]
dim(d.sub)

#because this dataset has such high sequencing, try cutting by mean abundance, instead
d.sub <- d.sub[, apply(d.sub, 2, mean) > mean.cutoff]
dim(d.sub)

#duplicate d.sub to analyze again at the end of script
h <- as.matrix(d.sub)

#First difference the data, since it's on a daily scale 
d.sub.first.diff <- d.sub[-dim(d.sub)[1], ] - d.sub[-1, ]
dim(d.sub.first.diff) #same number of OTUs
  
cor.mat.true <- cor(d.sub.first.diff)

par(mfrow = c(1,1))

d.sub <- as.matrix(d.sub)

otu.cors.true.pos <- apply(cor.mat.true, 1, pos.mean)
otu.cors.true.neg <- apply(cor.mat.true, 1, neg.mean)

#define number of permutations for null model 
perm.cor <- 50 #effect sizes for 200 are within 0.01

#create vector to hold median otu-otu correlations for each otu
med.otu.cors <- vector()

#create random number vector to hold seeds for randomization
#set seed here to create same seeds in each run
set.seed(888)
seeds <- sample(seq(1:100000000), size = (dim(d.sub)[2] * perm.cor * dim(d.sub)[1]))

if(!col.shuffle){
  #This is the row-shuffle null model
  for(m in 1:dim(d.sub)[2]){
    which.otu <- m
    
    #create vector to hold correlations from every permutation for each single otu
    perm.cor.vec.mat <- vector()
    
    #Run this loop as many times as specified by perm.cor
    for(i in 1:perm.cor){
      #Duplicate the d.sub matrix, and then replace entries with randomized numbers
      d.sub2 <- d.sub
      for(j in 1:dim(d.sub)[1]){ 
        #Randomize only taxa present in sample (i.e. abundance > 0)
        which.replace <- which(d.sub2[j, ] > 0 ) 
        #Do not randomize abundance of focal taxon
        if(d.sub2[j, m] > 0) {which.replace <- which.replace[!which.replace == m]}
        which.replace.minus <- which.replace[!(which.replace %in% m)]
        #Set new seed for upcoming randomization
        set.seed(seeds[j + (i-1)*(dim(d.sub)[1]) + (m-1)*dim(d.sub)[2] ])
        #Replace original values with randomized values
        d.sub2[j, which.replace.minus ] <- sample(d.sub[ j, which.replace.minus]) 
        
      }
      #replace focal column with original column, just to be sure it stays the same
      d.sub2[, which.otu] <- d.sub[ , which.otu]
      
      #calcuinv correlations between randomized taxon vectors
      d.sub2.first.diff <- d.sub2[-dim(d.sub2)[1], ] - d.sub2[-1, ]
      cor.mat.true.null <- cor(d.sub2.first.diff)
      
      #save the vector corresponding to the focal taxon, m
      perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.true.null[,m])
      
    }
    #Take the median of the correlations produced by the null model
    med.otu.cors <- cbind(med.otu.cors, apply(perm.cor.vec.mat, 1, median))
    
    if(m %% 20 == 0){print(m)}
  }
  
} else {
  #This is the column shuffle null model
  for(m in 1:dim(d.sub)[2]){
    which.otu <- m
    
    #create vector to hold correlations from every permutation for each single otu
    perm.cor.vec.mat <- vector()
    
    #Run this loop as many times as specified by perm.cor
    for(i in 1:perm.cor){
      #Duplicate the d.sub matrix, and then replace entries with randomized numbers
      d.sub2 <- d.sub
      for(j in 1:dim(d.sub)[2]){ 
        #Set new seed for upcoming randomization
        set.seed(seeds[j + (i-1)*(dim(d.sub)[1]) + (m-1)*dim(d.sub)[2] ])
        
        #randomize each taxon's abundance vector
        d.sub2[, j ] <- sample(d.sub[ ,j ]) 
      }
      #replace focal column with original column
      d.sub2[, which.otu] <- d.sub[ , which.otu]
      
      #calcuinv correlations between randomized taxon vectors
      d.sub2.first.diff <- d.sub2[-dim(d.sub2)[1], ] - d.sub2[-1, ]
      cor.mat.true.null <- cor(d.sub2.first.diff)
      
      #save the vector corresponding to the focal taxon, m
      perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.true.null[,m])
      
    }
    #Take the median of the correlations produced by the null model
    med.otu.cors <- cbind(med.otu.cors, apply(perm.cor.vec.mat, 1, median))
    
    if(m %% 20 == 0){print(m)}
  }
}

#get observed - expected individual correlations
obs.exp.cors.each <- cor.mat.true - med.otu.cors
diag(obs.exp.cors.each) <- 0

#calcuinv connectedness values
obs.exp.true.pos <- apply(obs.exp.cors.each, 2, pos.mean)
obs.exp.true.neg  <- apply(obs.exp.cors.each, 2, neg.mean) 

#####################################################################################################
#####################################################################################################

#duplicate full matrix (after detection limit imposed) to do further analyses 
a.full <- b

#create full pairwise bray curtis dissimilarity matrix 
bc <- as.matrix(vegdist(a.full))

# #pull out turnover from time t to t+1
# bc.diff <- vector()
# for(i in 1:(dim(h)[1]-1)){
#   bc.diff[i] <- bc[i, i+1]
# }
# bc.diff.all <- bc.diff

coh.neg <- d.sub %*% obs.exp.true.neg
coh.pos <- d.sub %*% obs.exp.true.pos

#####################################################################################################
#####################################################################################################

#Run in a loop to generate fig with null model 
par(mfrow = c(2, 1))
for(TF in c(F, T)) { 
  
  ######################### Define params for invasion analysis: 
  
  # Find the number of taxa that are absent until at least day n of the time series
  streak.cutoff <- 0.01
  tax.pers.cutoff <- .5
  use.rare.noninv <- TF #if T, uses uncommon but non-invasive taxa
  nonstreak.cutoff <- .3
  
  q <- a.full[, colSums(a.full) > 0] #take out zero columns
  q <- q[, apply(ceiling(q), 2, sum) <= dim(q)[1] * tax.pers.cutoff]
  dim(q)
  
  streak.pvals <- as.numeric(apply(q, 2, streakp))
  max.streak <- apply(q, 2, streakval)
  
  #For any taxa present in all samples, set p value to 1
  streak.pvals[apply(ceiling(q), 2, mean) == 1] <- 1
  
  tax.pers <- apply(ceiling(q), 2, mean)
  #hist(max.streak[tax.pers < tax.pers.cutoff])
  
  #Identify invasive taxa 
  inv <- streak.pvals < streak.cutoff  & tax.pers < tax.pers.cutoff
  if(use.rare.noninv){ inv <- streak.pvals > nonstreak.cutoff  & tax.pers < tax.pers.cutoff }
  sum(inv)
  mean(rowSums(q[, inv]))
  
  #create presence absence matrix to see how long these taxa persist
  q.pa <- ceiling(q)
  new.pers <- colSums(q.pa[, inv])
  #hist(new.pers, breaks = 20)
  #how many taxa persist for at least 10 days?
  sum(new.pers > 10) 
  #how many taxa only show up once?
  sum(new.pers == 1) 
  #make a histogram of the max abundance attained by these taxa
  #hist(log(apply(q[, inv], 2, max)), breaks = 200)
  
  inv.sum <- apply(q[, inv], 1, sum) 
  # plot(inv.sum)
  # log.inv <- log(rowSums(q[, inv]))
  # plot(log.inv)
  
  ####################################
  
  #Order by most abundant "invader" taxa
  q.inv <- q[, inv]
  q.ord <- q.inv[, order(colSums(q.inv), decreasing = T)]
  
  #Take out taxa where there is only 1 observation
  q.inv <- q.inv[, apply(ceiling(q.inv), 2, sum) > 1]
  dim(q.inv)
  mean(rowSums(q.inv))
  
  ## Make a model of each invader taxon's abundance as a fxn of previous abundance and cohesion
  
  w.inv <- q.inv 
  which.inv.t <- which(w.inv[-dim(w.inv)[1], ] > 0, arr.ind = T)
  which.inv.t1 <- which.inv.t
  which.inv.t1[, 1] <- which.inv.t[, 1] + 1 #should be + 1
  
  inv.vec.t <- w.inv[which.inv.t]
  inv.vec.t1 <- w.inv[which.inv.t1] 
  otu.id <- as.factor(which.inv.t[, 2])
  
  c.neg <- d.sub %*% obs.exp.true.neg
  c.pos <- d.sub %*% obs.exp.true.pos
  
  c.neg <- c.neg[rownames(d.sub) %in% rownames(w.inv)]
  c.pos <- c.pos[rownames(d.sub) %in% rownames(w.inv)]
  
  c.neg.inv <- c.neg[which.inv.t[,1]]
  c.pos.inv <- c.pos[which.inv.t[,1]]
  
  ########################################################################
  ########################################################################
  
  #What about subsetting to just instances where taxa are present for two consecutive times?
  consec <- inv.vec.t1 > 0
  
  inv.vec.con.t1 <- inv.vec.t1[consec]
  inv.vec.con.t <- inv.vec.t[consec]
  c.neg.inv.con <- c.neg.inv[consec]
  c.pos.inv.con <- c.pos.inv[consec]
  otu.id.con <- otu.id[consec]
  
  #Create competing models to evaluate fit via change in AIC
  lm.cons.1 <- (lm(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con))
  lmer.cons.1 <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (1|otu.id.con))
  lmer.cons.1.1 <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (1 + c.pos.inv.con   |otu.id.con), REML = F)
  lmer.cons.1.2 <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (1 + log(inv.vec.con.t) + c.pos.inv.con   |otu.id.con), REML = F)
  lmer.cons.1.2a <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (0 + log(inv.vec.con.t) + c.pos.inv.con |otu.id.con), REML = F)
  lmer.cons.1.2b <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (1 +  c.pos.inv.con   |otu.id.con), REML = F)
  lmer.cons.1.2c <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (1 + log(inv.vec.con.t)   |otu.id.con), REML = F)
  log.inv.vec.con.t <- log(inv.vec.con.t)
  
  summary(lm.cons.1) 
  summary(lmer.cons.1) 
  summary(lmer.cons.1.1)
  summary(lmer.cons.1.2)
  
  AIC(lm.cons.1)
  AIC(lmer.cons.1) 
  AIC(lmer.cons.1.1)
  AIC(lmer.cons.1.2)
  AIC(lmer.cons.1.2) - AIC(lmer.cons.1.2a)
  AIC(lmer.cons.1.2) - AIC(lmer.cons.1.2b)
  AIC(lmer.cons.1.2) - AIC(lmer.cons.1.2c)
  AIC(lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t)  + (1|otu.id.con))) # cohesion improves model fit 
  
  r.squaredGLMM(lmer.cons.1.1) 
  r.squaredGLMM(lmer.cons.1.2)
  
  #find median peak abundance
  median(apply(q[, inv], 2, max))
  #find IQR peak abundance 
  sort(apply(q[, inv], 2, max))[c(round(sum(inv) * .25), round(sum(inv) * .75))]
  #find median streak length
  median(max.streak[inv])
  #Get IQR of strek lengths 
  sort(max.streak[inv])[c(round(sum(inv) * .25), round(sum(inv) * .75))]
  #Find average proportion of community that is invader
  mean(rowSums(q.inv))
  sum(inv)
  
  
  otu.id.factor <- as.factor(otu.id.con)
  
  #Plot fitted lines for relatively abundant invasive and non-invasive taxa
  otu.abun <- .01
  otu.seq <- 50
  #par(mfrow = c(1, 1))
  #Remake MLM plot using real estimated slopes/intercepts from n most abundant taxa 
  how.many.abun <- 50
  q.inv.model <- q.inv[, unique(otu.id.con)]
  dim(q.inv.model)
  most.ab <- which(log(apply(q.inv, 2, mean)) %in% sort(log(apply(q.inv.model, 2, mean)), decreasing = T)[1:how.many.abun] ) 
  plot(1, 1, col = "white", xlim = range(c.pos.inv), ylim = c(.0001, .1), xlab = "Positive Cohesion", ylab = "Predicted Abundance", log = "y", bty = "l")
  for(ab.otu in most.ab){
    coh.min <- min(c.pos.inv.con[as.numeric(otu.id.con) == ab.otu])
    coh.max <- max(c.pos.inv.con[as.numeric(otu.id.con) == ab.otu])
    otu.int <- ranef(lmer.cons.1.1)$otu.id.con[, 1][which(unique(otu.id.con) == ab.otu )] + fixef(lmer.cons.1.1)[1]
    otu.slope <- ranef(lmer.cons.1.1)$otu.id.con[, 2][which(unique(otu.id.con) == ab.otu )]
    otu.y0 <- otu.int + otu.slope*coh.min
    otu.y1 <- otu.int + otu.slope*coh.max
    
    otu.newdata <- data.frame( inv.vec.con.t = rep(otu.abun, times = 2), c.pos.inv.con = c(coh.min, coh.max), c.neg.inv.con = rep(-0.08, times = 2), otu.id.con = rep(ab.otu, times = 2) )
    
      pred.otu <- predict(lmer.cons.1.2, newdata = otu.newdata, re.form = NULL)
    
    if(exp(pred.otu[2]) > 1.5 * exp(pred.otu[1] )){
      plot.col <- "blue"
    } else{
      if(exp(pred.otu[2]) * 1.5 < exp(pred.otu[1] ) ){
        plot.col <- "red"
      } else{
        plot.col <- "black"
      }
    }
    segments(x0 = coh.min, x1 = coh.max, y0 = exp(pred.otu[1]), y1 = exp(pred.otu[2]), col = plot.col, lwd = .8)  
  }
  
  print(summary(lmer.cons.1.1))
  print(summary(lmer.cons.1.2))
  
}

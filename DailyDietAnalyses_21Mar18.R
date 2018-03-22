## Analysis to identify and study invasive OTUs in the human gut
## 30 Jan 18 

#library(reshape)
#library(Hmisc)
#library(plyr)
#library(vegan)
#library(AICcmodavg)
#library(devtools)
#library(car)
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

trim <- read.csv("daily_diet_trim_meta_25Jul17.csv", header = T)
meta <- read.csv("DailyDietMetaStool_25Jul17.csv", header = T)

#take out metadata, which is the first 2 rows
ab <- trim[-c(1,2), -1]
rownames(ab) <- paste(rep("otu", times = dim(ab)[1]), seq(1, dim(ab)[1], 1), sep = ".")
ab <- t(ab)

############################### input options for analysis
#pick a minimum proportion of samples to analyze OTU
zero.prop <- .5 #usually .5
#choose mean abundance cutoff
mean.cutoff <- 0.001 
col.shuffle <- F #(T also gives similar results)
subject <- "A"
det.cutoff <- -10
use.detection.limit <- T
randomize.rows <- F #to randomize values before calculating cohesion for null model

#########################################################
#take out stool 69 since it has no metadata
ab <- ab[rownames(ab) != "Stool69", ]

##create vectors of host ID and day ID
days.all <- meta$COLDAY
host.all <- unlist(lapply(meta$HOST, substr, start = 11, stop = 11))
#put the vectors in the order of the samples 
days <- days.all[match(rownames(ab), meta$ANONYMIZED_NAME)]
host <- host.all[match(rownames(ab), meta$ANONYMIZED_NAME)]
length(days)
dim(ab)

#Subset to specified host and correct day range.
if(subject == "B") { a <- ab[host == "B" & days < 250 , ] 
day.ord <- days[host == "B" & days < 250 ]} #sampling is very sparse after day 250 
if(subject == "A") { a <- ab[host == "A"  , ] 
day.ord <- days[host == "A"  ]} 

#Put into correct order by day
a <- a[order(day.ord), ]

colnames.save <- colnames(a)
rownames.save <- rownames(a)

#put into matrix
a.mat <- matrix(as.numeric(a), nrow = nrow(a), ncol = ncol(a))
a <- a.mat
colnames(a) <- colnames.save
rownames(a) <- rownames.save
dim(a)

#Now, create global detection limit
if(use.detection.limit) {
det.mat <- a
det.mat[det.mat == 0] <- 1
hist(log(apply(det.mat, 1, min)), breaks = 50)

#Take out samples with detection lower than universal detection cutoff
a <- a[log(apply(det.mat, 1, min)) < det.cutoff, ]
a[a < exp(det.cutoff)] <- 0
hist(rowSums(a))
}

#b <- a

#create new matrix to manpiulate
d.sub <- as.matrix(a)

#take out OTUs with fewer occurrences than the presence cutoff 
## This has usually been commented out for this anlaysis. Results very similar if enabled
zero.cutoff <- ceiling(zero.prop*dim(d.sub)[1])
d.sub <- d.sub[, apply(d.sub, 2, zero) < (dim(d.sub)[1]-zero.cutoff) ]
#dim(d.sub)

#because this dataset has such high sequencing, cut by mean abundance, also
d.sub <- d.sub[, apply(d.sub, 2, mean) > mean.cutoff]
dim(d.sub)


# For randomizing values before calculating cohesion
if(randomize.rows) {
  for(row in 1:dim(d.sub)[1]){
   d.sub[row, ] <- sample(d.sub[row, ])
  }
}

#First difference the data, since it's on a daily scale 
d.sub.first.diff <- d.sub[-dim(d.sub)[1], ] - d.sub[-1, ]
dim(d.sub.first.diff) #same number of OTUs

#Look at difference in autocorrelation after first differencing
# cor.after.diff <- vector()
# for(tx in 1:dim(d.sub.first.diff)[2]){
# cor.after.diff[tx]  <- cor(d.sub.first.diff[ -dim(d.sub.first.diff)[1] , tx] , d.sub.first.diff[ -1 , tx ])
# }
# hist(cor.after.diff)
# 
# cor.before.diff <- vector()
# for(tx in 1:dim(d.sub)[2]){
#   cor.before.diff[tx]  <- cor(d.sub[ -dim(d.sub)[1] , tx] , d.sub[ -1 , tx ])
# }
# hist(cor.before.diff)
  
cor.mat.true <- cor(d.sub.first.diff)

par(mfrow = c(1,1))

d.sub <- as.matrix(d.sub)

otu.cors.true.pos <- apply(cor.mat.true, 1, pos.mean)
otu.cors.true.neg <- apply(cor.mat.true, 1, neg.mean)

#define number of permutations for null model 
perm.cor <- 200

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
      
      #calculate correlations between randomized taxon vectors
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
      
      #calculate correlations between randomized taxon vectors
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

#calculate connectedness values
obs.exp.true.pos <- apply(obs.exp.cors.each, 2, pos.mean)
obs.exp.true.neg  <- apply(obs.exp.cors.each, 2, neg.mean) 

#####################################################################################################
#####################################################################################################

#duplicate full matrix (after detection limit imposed) to do further analyses 
a.full <- a 

coh.neg <- d.sub %*% obs.exp.true.neg
coh.pos <- d.sub %*% obs.exp.true.pos

#####################################################################################################
#####################################################################################################

#Run in a loop to generate figure of mlm results of invasive and non-invasive OTUs vs positive cohesion
par(mfrow = c(2, 1))
for(TF in c(F, T)) { 

######################### Define params for invasion analysis: 

# Find the number of taxa that are absent until at least day n of the time series
streak.cutoff <- 0.01
tax.pers.cutoff <- .5 #.3 shows exact same pattern for A
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

#rownames(q.inv)
w.inv <- q.inv 
#Find indices of all instances where a taxon is present. 
which.inv.t <- which(w.inv[-dim(w.inv)[1], ] > 0, arr.ind = T)
which.inv.t1 <- which.inv.t
#Find indices of all instances immediately AFTER instances where the taxa are present
which.inv.t1[, 1] <- which.inv.t[, 1] + 1 #should be + 1

#Pull out non-zero abundances
inv.vec.t <- w.inv[which.inv.t]
#Double check all abundances are non-zero
sum(inv.vec.t == 0)
#Pull out abundances at next time point
inv.vec.t1 <- w.inv[which.inv.t1] 
#Find OTU identifier
otu.id <- as.factor(which.inv.t[, 2])

c.neg <- d.sub %*% obs.exp.true.neg
c.pos <- d.sub %*% obs.exp.true.pos

#removes cohesion values for any samples where there are no invasive taxa
#c.neg <- c.neg[rownames(d.sub) %in% rownames(w.inv)]
#c.pos <- c.pos[rownames(d.sub) %in% rownames(w.inv)]

c.neg.inv <- c.neg[which.inv.t[,1]]
c.pos.inv <- c.pos[which.inv.t[,1]]

#Check that vectors for cohesion and taxon abundances are the same length
length(c.neg.inv) == length(inv.vec.t)


########################################################################
########################################################################

#Find the instances where the abundance of the taxon at the consecutive time point is greater than 0
consec <- inv.vec.t1 > 0

#Create variables corresponding to only time points with consecutive presences
inv.vec.con.t1 <- inv.vec.t1[consec]
inv.vec.con.t <- inv.vec.t[consec]
c.neg.inv.con <- c.neg.inv[consec]
c.pos.inv.con <- c.pos.inv[consec]
otu.id.con <- otu.id[consec]

#Create competing models to evaluate fit via change in AIC
lm.cons.1 <- (lm(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con))
lmer.cons.1 <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (1|otu.id.con))
lmer.cons.1.1 <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (1 + c.pos.inv.con   |otu.id.con), REML = F)
lmer.cons.1.1a <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (0 + c.pos.inv.con   |otu.id.con), REML = F)
lmer.cons.1.1b <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (1  |otu.id.con), REML = F)
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
AIC(lmer.cons.1.1) - AIC(lmer.cons.1.1a)
AIC(lmer.cons.1.1) - AIC(lmer.cons.1.1b)
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

##Plot fitted lines for relatively abundant invasive and non-invasive taxa
#Set initial abundance for the predictions
otu.abun <- .01

#Make MLM plot using real estimated slopes/intercepts from n most abundant taxa 
how.many.abun <- 50
#create matrix of the OTUs modeled in the mlm. this removed non-invasive OTUs that were not present in consecutive days, and therefore were not included in the mlm
q.inv.model <- q.inv[, unique(otu.id.con)]
dim(q.inv.model)
#find the OTUs that are most abundant among the modeled OTUs
most.ab <- which(log(apply(q.inv, 2, mean)) %in% sort(log(apply(q.inv.model, 2, mean)), decreasing = T)[1:how.many.abun] ) 
#Make an empty plot
plot(1, 1, col = "white", xlim = range(c.pos.inv), ylim = c(.0001, .1), xlab = "Positive Cohesion", ylab = "Predicted Abundance", log = "y", bty = "l")
#for each OTU in the matrix of abundant OTUs, plot the fitted abundance against cohesion
for(ab.otu in most.ab){
  coh.min <- min(c.pos.inv.con[as.numeric(otu.id.con) == ab.otu])
  coh.max <- max(c.pos.inv.con[as.numeric(otu.id.con) == ab.otu])
  otu.int <- ranef(lmer.cons.1.1)$otu.id.con[, 1][which(unique(otu.id.con) == ab.otu )] + fixef(lmer.cons.1.1)[1]
  otu.slope <- ranef(lmer.cons.1.1)$otu.id.con[, 2][which(unique(otu.id.con) == ab.otu )]
  otu.y0 <- otu.int + otu.slope*coh.min
  otu.y1 <- otu.int + otu.slope*coh.max
  
  otu.newdata <- data.frame( inv.vec.con.t = rep(otu.abun, times = 2), c.pos.inv.con = c(coh.min, coh.max), c.neg.inv.con = rep(-0.08, times = 2), otu.id.con = rep(ab.otu, times = 2) )
  
  if(TF == T & subject == "A"){
     pred.otu <- predict(lmer.cons.1.1, newdata = otu.newdata, re.form = NULL)
  } else {
    pred.otu <- predict(lmer.cons.1.2, newdata = otu.newdata, re.form = NULL)
  }
  
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

#calculate change in expected abundance based on a 0.05 decrease in positive cohesion
#exp(1 + .05 * 6.28) / exp(1)
#exp(1 + .05 * 6.45) / exp(1)
#exp(1 + .05 * 8.66) / exp(1)

}


if(subject == "A" & TF == F){
  #Make 4 plots of invaders and non-invaders
  noninv <- streak.pvals > nonstreak.cutoff
  q.noninv <- q[, noninv]
  dim(q.noninv)
  q.noninv <- q.noninv[, order(colSums(q.noninv), decreasing = T)]
  
  par(mfrow = c(4, 1))
  colnames(q.noninv)[28]
  plot(q.noninv[, 14], xlab = "Day", ylab = "Relative Abundance", main = "", cex = .65, cex.axis = 1.3, cex.lab = 1.3, ylim = c(0, 1.4* max(q.noninv[, 14])), pch = 16, bty = "l", xaxs="i" )
  #text(310, 1.2*max(q.noninv[, 28]), "a", cex = 1.4)
  text(paste0("OTU 1: p = ", round(streak.pvals[match(colnames(q.noninv)[14], colnames(q))], 2) ) , x = 160, y = 1.3 *max(q.noninv[, 14]), cex = 1.2)
  colnames(q.noninv)[47]
  plot(q.noninv[, 54] , xlab = "Day", ylab = "Relative Abundance", main = "", cex = .65, cex.axis = 1.3, cex.lab = 1.3, ylim = c(0, 1.3* max(q.noninv[, 54])), pch = 16, bty = "l", xaxs="i" )
  #text(310, 1.1*max(q.noninv[, 57]), "b", cex = 1.4)
  text(paste0("OTU 2: p = ", round(streak.pvals[match(colnames(q.noninv)[54], colnames(q))], 2) ) , x = 160, y = 1.2 *  max(q.noninv[, 54]), cex = 1.2)
  colnames(q.ord)[2]
  plot(q.ord[, 2] , xlab = "Day", ylab = "Relative Abundance", main = "", cex = .65, cex.axis = 1.3, cex.lab = 1.3, ylim = c(0, 1.3* max(q.ord[, 2])), pch = 16, bty = "l", xaxs="i" )
  #text(310, 1.1*max(q.ord[, 2]), "c", cex = 1.4)
  text("OTU 3: p < 0.001" , x = 165, y = 1.2 * max(q.ord[, 2]), cex = 1.2)
  colnames(q.ord)[26]
  plot(q.ord[, 26] , xlab = "Day", ylab = "Relative Abundance", main = "", cex = .65, cex.axis = 1.3, cex.lab = 1.3, ylim = c(0, 1.3* max(q.ord[, 26])) , pch = 16, bty = "l" , xaxs="i" )
  #text(310, 1.1*max(q.ord[, 26]), "d", cex = 1.4)
  text("OTU 4: p < 0.001" , x = 165, y =  1.2 * max(q.ord[, 26]), cex = 1.2)
}


#As additional null model, randomize connectedness values and see how R2 values compare to true 
if(TF == F){
  r2.fixef <- vector()
  pos.coh.coef <- vector()
  for(it in 1:1000){
    shuffle.ord <- sample(seq(1:dim(d.sub)[2]))
    c.neg <- d.sub %*% obs.exp.true.neg[shuffle.ord]
    c.pos <- d.sub %*% obs.exp.true.pos[shuffle.ord]
    
    #c.neg <- d.sub %*% sample(obs.exp.true.neg, replace = T)
    #c.pos <- d.sub %*% sample(obs.exp.true.pos, replace = T)
    
    c.neg <- c.neg[rownames(d.sub) %in% rownames(w.inv)]
    c.pos <- c.pos[rownames(d.sub) %in% rownames(w.inv)]
    
    c.neg.inv <- c.neg[which.inv.t[,1]]
    c.pos.inv <- c.pos[which.inv.t[,1]]
    
    c.neg.inv.con <- c.neg.inv[consec]
    c.pos.inv.con <- c.pos.inv[consec]
    
    #plot(c.pos ~ d.sub %*% obs.exp.true.pos)
    
    lmer.cons.iter <- lmer(log(inv.vec.con.t1) ~ log(inv.vec.con.t) + c.neg.inv.con + c.pos.inv.con + (1  + c.pos.inv.con   |otu.id.con), REML = F)
    #summary(lmer.cons.iter)
    r2.fixef[it] <- r.squaredGLMM(lmer.cons.iter)[1]
    pos.coh.coef[it] <- summary(lmer.cons.iter)$coefficients[4, 1] / summary(lmer.cons.iter)$coefficients[4, 2]

    print(it)
  }
  par(mfrow = c(1, 1))
  hist(r2.fixef)
  hist(pos.coh.coef)
  mean(abs(pos.coh.coef) > 6.16 ) 
}




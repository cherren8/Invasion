#!/usr/bin/env Rscript

#Arguments input should be competition coefficient, exchange proportion, and number of runs

args = commandArgs(trailingOnly=TRUE)

## 06 Nov 17
## Modeling invasion using metabolite explicit model 
library(combinat)
library(data.table)
#library(vegan)
library(abind)
#library(ggplot2)
library(dplyr)
library(scales)

## Create necessary functions
min.of.n <- function(vec, min.n){
  sort.vec <- sort(vec, decreasing = T)
  min.sub <- min(sort.vec[1:min.n])
  return(min.sub)
}



#Variable inputs based on command line arguments
comp.mean <- as.numeric(args[1])
#print(comp.mean)
comp.sd <- comp.mean * .3    
direct.prop <- as.numeric(args[2])
num.runs <- as.numeric(args[3])

result.list <- list()

## Fixed params
min.comp <- .01
inv.comp <- .9
input.rate <- 100
num.taxa <- 15 
num.all.taxa <- num.taxa + 1
num.mets.all <- 7 
mets.req <- 4
inv.mets.req <- 4 
mets.ex <- 2
flush <- .15
time.steps <- 20000 #if the model doesn't equilibrate in this time, it is restarted
input.all.mets <- F
flow.strength <- 1
over.step.count <- 0


#create a vector to record number of taxa coexisting
num.coex <- vector()
#create a vector to hold equilibrium met pools
eq.mets <- vector()
#Create a vector for sum of metabolites
mets.all <- vector()
#create a vector to record number of taxa coexisting BEFORE invader
num.coex.preinv <- vector()
#create a vector to hold equilibrium met pools BEFORE invader
eq.mets.preinv <- vector()
#create a vector to look at whether invader persists
inv.persist <- vector()
#make a vector to hold final abundances of taxa before invasion
final.abuns.preinv <- vector()
#make a vector to hold final abundances of taxa 
final.abuns <- vector()

#############################################################################################
#############################################################################################

run <- 1
while(run <= num.runs ){
  
  reqs <- vector()
  n.mets.for.taxa <- vector()
  
  for(k in 1:length(mets.req)){
    all.mets.k <- unique(permn(c(rep(1, mets.req[k]), c(rep(0, num.mets.all - mets.req[k])))))
    which.mets.k <- sample(seq(1, length(all.mets.k)), round(num.taxa / length(mets.req) ), replace = F)
    
    reqs.k <- matrix(unlist(all.mets.k[which.mets.k]), ncol = round(num.taxa / length(mets.req) ), byrow = F)
    
    reqs <- cbind(reqs, reqs.k)
    dim(reqs)  
    
    n.mets.for.taxa <- c(n.mets.for.taxa, rep(mets.req[k], round(num.taxa / length(mets.req))) )
  }
  n.mets.for.taxa <- c(n.mets.for.taxa, inv.mets.req)
  
  #define invader reqs
  all.inv.reqs <- unique(permn(c(rep(1, inv.mets.req), c(rep(0, num.mets.all - inv.mets.req)))))
  all.inv.reqs.mat <- matrix(unlist(all.inv.reqs), nrow = num.mets.all , byrow = F)
  M1 = setkey(data.table(t(reqs)))
  
  M2 = setkey(data.table(t(all.inv.reqs.mat)))
  shared <- na.omit(
    M2[M1,which=TRUE]
  )
  shared <- as.vector(shared)
  ifelse(length(shared) >= 1, 
         which.inv.reqs <- sample(seq(1, dim(M2)[1], 1)[-shared], 1), 
         which.inv.reqs <- sample(seq(1, dim(M2)[1], 1), 1) )
  inv.reqs <- M2[which.inv.reqs, ]
  
  reqs <- cbind(reqs, t(inv.reqs))
  
  dim(reqs)
  
  
  #mets.in <- unlist(sample(unique(permn(c(rep(1, max(mets.req) ), c(rep(0, num.mets.all - max(mets.req)))))), 1))
  mets.in <- reqs[, sample(seq(1, dim(reqs)[2] - 1, 1), 1)]
  if(input.all.mets) mets.in <- rep(1, times = num.mets.all)
  
  all.ex.poss <- matrix(as.numeric(!reqs), ncol = num.all.taxa)
  how.many.ones <- function(vec, ones){
    which.ones <- which(vec == 1)
    keep.ones <- sample(which.ones, ones)
    vec.zero <- rep(0, length(vec))
    vec.zero[keep.ones] <- 1
    return(vec.zero)
  }
  ex <- apply(all.ex.poss, 2, how.many.ones, ones = mets.ex)
  
  #Make a matrix holding the donor, recipient, and met for each possible directed flow 
  donors <- vector()
  recipients <- vector()
  met.flow <- vector()
  for(q in 1:num.mets.all){
    for(w in 1:num.all.taxa){
      if(ex[q, w] == 1){
        num.req <- sum(reqs[q, ])
        donors <- c(donors, rep(w, num.req))
        recipients <- c(recipients, which(reqs[q, ] == 1))
        met.flow <- c(met.flow, rep(q, times = num.req))
      }
    }
  }
  possible.flows <- cbind(donors, recipients, met.flow)
  possible.flows <- possible.flows[possible.flows[,2] != num.all.taxa, ]
  possible.flows <- possible.flows[possible.flows[,1] != num.all.taxa, ]
  
  real.flows <- possible.flows[sample(seq(1, dim(possible.flows)[1]), ceiling(direct.prop*dim(possible.flows)[1])), ]
  real.flows <- as.matrix(real.flows)
  dim(real.flows)
  
  #Put flows into an array of num.taxa x num.taxa x num.mets 
  flow.array <- array(0, dim = c(num.all.taxa, num.all.taxa, num.mets.all))
  donor.array <- array(0, dim = c(num.all.taxa, num.all.taxa, num.mets.all))
  for(which.met in 1:num.mets.all){
    donor.array[real.flows[, 1][real.flows[, 3] == which.met], real.flows[, 2][real.flows[, 3] == which.met], which.met] <- -1
    flow.row <- real.flows[, 2][real.flows[, 3] == which.met]
    flow.col <- real.flows[, 1][real.flows[, 3] == which.met]
    for(pair in 1:length(flow.row)){
      flow.array[flow.row[pair], flow.col[pair], which.met] <- 1
    }
  }
  
  ## Variables 
  abun.series <- vector()
  mets.series <- vector()
  
  ##  Initialize vectors
  init.abun <- 50
  abuns <- c(rep(init.abun, num.taxa ), 0)
  mets.pool <- rep(0, num.mets.all)
  met.stores <- matrix(0, num.mets.all, num.all.taxa)
  
  #pull competitors from a uniform distribution but give invader average competitive value
  comp.vec <- c(rnorm(num.taxa, comp.mean, comp.sd) , inv.comp)   #rnorm(num.all.taxa, comp.mean, comp.sd)
  comp.vec <- apply(cbind(comp.vec, rep(min.comp, times = length(comp.vec))), 1, max)
  comp <- matrix(rep(comp.vec, num.mets.all), ncol = num.all.taxa, byrow = T)
  dim(comp)
  
  stop.val <- 0.01
  max.abun.change <- 1
  i <- 1
  inv.added <- F
  
  while(!(max.abun.change <= stop.val | i >= time.steps)){ #if either of these are true, loop stops 
    mets.pool <- mets.pool  +  (input.rate * mets.in)
    
    demand.mat <-  (reqs * comp)  *  ( matrix(rep(abuns, each = num.mets.all), num.mets.all, num.all.taxa) - met.stores ) 
    demand <-  rowSums ( (reqs * comp)  *  ( matrix(rep(abuns, each = num.mets.all), num.mets.all, num.all.taxa) - met.stores ) )
    
    comp.abil <- (matrix(rep(abuns, each = num.mets.all), num.mets.all, num.all.taxa) - met.stores ) * reqs * comp
    
    comp.abil.rel <- comp.abil / rowSums(comp.abil)
    comp.abil.rel[is.nan(comp.abil.rel)] <- 0 #for cases where row sums are 0, will return infinity. should be 0
    
    mets.add <- (apply(cbind(demand, mets.pool), 1, min)) * comp.abil.rel
    
    met.stores <- met.stores  +  mets.add
    mets.pool <- mets.pool - rowSums(mets.add)
    mets.pool <- apply(cbind(rep(0, length(mets.pool)), mets.pool), 1, max)
    
    #Find the amount of only metabolites that are required by each taxon
    #req.met.stores <- matrix(met.stores[as.logical(reqs)] , max(mets.req), num.taxa)
    
    #growth <- apply(req.met.stores, 2, min)
    #growth[is.nan(growth)] <- 0
    
    #Make loop to calculate min of requried mets
    growth <- vector()
    for(nt in 1:num.all.taxa){
      growth.n <- min.of.n(met.stores[, nt], min.n = n.mets.for.taxa[nt])
      growth <- c(growth, growth.n)
    }
    
    abuns <- abuns + growth 
    
    met.ex <- t( growth * t(ex))
    
    mets.traded <- vector()
    for(met.i in 1:num.mets.all){
      flow.norm <- sweep(flow.array[,, met.i], 2, colSums(flow.array[,, met.i]), '/')
      flow.norm[is.na(flow.norm)] <- 0
      flow.norm <- flow.norm * flow.strength
      met.avail <- flow.norm %*% met.ex[met.i, ]
      met.missing <- abuns - met.stores[met.i, ]
      mets.transfer <- apply(cbind(met.avail, met.missing), 1, min)
      met.stores[met.i, ] <- met.stores[met.i, ]  +  mets.transfer
      mets.traded[met.i] <- sum(mets.transfer)
    }
    
    #Subtract exchanged mets from mets to be excreted
    met.stores <- (met.stores - t( growth * t(reqs)) ) * (1 - flush)
    mets.pool <- (mets.pool  +  rowSums(met.ex) - mets.traded) * (1 - flush)
    mets.pool <- apply(cbind(rep(0, length(mets.pool)), mets.pool), 1, max)
    
    #abuns <- apply(cbind(floor(abuns * (1 - flush)), rep(0, times = length(abuns))), 1, max) #makes taxa go extinct 
    abuns <- abuns * (1 - flush)
    
    abun.series <- cbind(abun.series, abuns)
    mets.series <- cbind(mets.series, mets.pool)
    
    #if(i %% 500 == 0) abuns[abuns > 1] <- sample(abuns[abuns > 1]) #perturb the system to see if it is stable
    
    if(i > 1) max.abun.change <- max(abs(abuns - abun.series[, i - 1]))
    
    #after initial equilibration, add an invader
    eq.scenario <- inv.added == F & max.abun.change < stop.val
    half.time <- i == time.steps / 2
    if(eq.scenario | half.time){
      num.coex.half <- sum(abuns > 1)
      eq.mets.half <-  mets.pool 
      abuns[length(abuns)] <- init.abun
      final.abuns.half <- abuns
      inv.added <- T
      max.abun.change <- mean(abs(abuns - abun.series[, i - 1])) #makes it so loop continues
      if(i == time.steps / 2) over.step.count <- over.step.count + 1
    }
    
    i <- i  +  1 
    
  }
  
  # if(i < time.steps) {
  num.coex.preinv <- c(num.coex.preinv, num.coex.half)
  eq.mets.preinv <- cbind(eq.mets.preinv, eq.mets.half) 
  final.abuns.preinv <- cbind(final.abuns.preinv, final.abuns.half)
  
  final.abuns <- cbind(final.abuns, abuns)
  num.coex <- c(num.coex, sum(abuns > 1))
  eq.mets <- cbind(eq.mets, mets.pool) 
  mets.all <- c(mets.all, sum(mets.pool))
  
  abun.pa <- abun.series[, i - 1]
  abun.pa[abun.pa >= 1] <- 1
  abun.pa[abun.pa < 1] <- 0
  inv.persist <- c(inv.persist, abun.pa[length(abun.pa)])
  
  run <- run + 1
  # } else {
  #   run <- run 
  # }
  
  # if(run %% 10 == 0 ) print(c(run, i) )
  #print(c(run, i, st) )
  #if(sum(abuns > 1) > num.mets.all & direct.prop == 0) stop() # if getting more coexistence than mets, check why
  
}



params <- c(input.rate, num.taxa, num.mets.all, mets.req, inv.mets.req, mets.ex, flush, time.steps, num.runs, comp.mean, comp.sd, flow.strength)
names(params) <- c("input.rate", "num.taxa", "num.mets.all", "mets.req", "inv.mets.req", 'mets.ex', "flush", "time.steps", "num.runs", "comp.mean", "comp.sd", "flow.strength")
params

over.step.prop <- over.step.count / num.runs

#Create files to record median metabolite pools, median diversity, median number of individuals, and proportion successful invaders for each proportion 
summary.results <- cbind(rep(comp.mean, times = num.runs), rep(direct.prop, times = num.runs),  num.coex.preinv, num.coex, colSums(final.abuns.preinv), colSums(final.abuns), colSums(eq.mets.preinv), colSums(eq.mets), inv.persist)
colnames(summary.results) <- c("comp.mean", "direct.prop", "num.coex.preinv", "num.coex.post", "total.ind.preinv", "total.ind.post", "total.mets.pre", "total.mets.post", "inv.persist")
rownames(summary.results) <- seq(1, dim(summary.results)[1], 1)

summary.stats <- c(apply(summary.results[, 1:8], 2, median), mean(summary.results[, 9]), over.step.prop)
names(summary.stats) <- c(colnames(summary.results), "over.step.prop")

# write.csv(t( summary.stats) , "MetModelOutput.csv")

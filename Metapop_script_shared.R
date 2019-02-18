#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################
#Metapopulation model for Stockholm archipelago data, which has 10 islands, and Elsewhere (other sites, all in the Baltic)
#Age is used as state, with with 5 stages: 1-year,2-year,3-year,4-year, adults.
#Adults are in turn separated into recently ringed and previously ringed
#Birds are only ringed as 1-year and adults - recently ringed.
#Thus only one group of individuals
#Dataset formatted as individual histories
#m-array format in analysis
#Last update 2019-02-18 /MKA

#Model 7
#Estimating survival for 2-4-yr-olds jointly
#Same p structure as in model 6b
#p is encounter probability

#Model 8 
#Same structure for phi and p as in model 7
#Including a distance function for psi
#p is encounter probability

#Model 9 
#Removing the Elsewhere states
#p is encounter probability

#Model 10
#Two different mu for p, one 1996-2013, one 2014-2017

#Model 11
#A different p at site 9 and 10 when resights occur


########################
#######################


########################
#######################
#It should be possible to run the code on rows 39 - ~160 and then start skip many rows and continue where
#model 10 starts

rm(list = ls())

library(jagsUI)
library(gdata)

# Get the data - encounter histories of birds captured 1995-2017
enc <- read.csv("Both1995-2017_DatasetCopy20190205.csv", sep=",") 
#Encounter data produced with EncounterHistoryScript2018.R
dist <- read.xls("ColonyProtectionSizeDistance.xlsx", sheet="Distance", na.strings=c(".", "", " ", "  ", "NA", "NULL", "__"), fileEncoding="latin1")
covar <- read.xls("ColonyProtectionSizeDistance.xlsx", sheet="Protection2018_Size", na.strings=c(".", "", " ", "  ", "NA", "NULL", "__"), fileEncoding="latin1")
covar$SiteNo. <- c(10,7,9,8,6,3,5,4,1,2) #The numbers used for the islands in the model structure

#Removing the birds ringed in 2017
enc <- enc[which(enc$Myear!=2017),]

# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
  # Unobservable: number of state that is unobservable
  n.occasions <- dim(PSI.STATE)[4] + 1
  CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
  g <- colSums(marked)
  for (s in 1:dim(PSI.STATE)[1]){
    if (g[s]==0) next  # To avoid error message if nothing to replace
    mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
  } #s
  for (i in 1:sum(marked)){
    for (s in 1:dim(PSI.STATE)[1]){
      if (mark.occ[i,s]==0) next
      first <- mark.occ[i,s]
      CH[i,first] <- s
      CH.TRUE[i,first] <- s
    } #s
    for (t in (first+1):n.occasions){
      # Multinomial trials for state transitions
      if (first==n.occasions) next
      state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
      CH.TRUE[i,t] <- state
      # Multinomial trials for observation process
      event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
      CH[i,t] <- event
    } #t
  } #i
  # Replace the NA and the highest state number (dead) in the file by 0
  CH[is.na(CH)] <- 0
  CH[CH==dim(PSI.STATE)[1]] <- 0
  CH[CH==unobservable] <- 0
  id <- numeric(0)
  for (i in 1:dim(CH)[1]){
    z <- min(which(CH[i,]!=0))
    ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
  }
  return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
  # CH: capture histories to be used
  # CH.TRUE: capture histories with perfect observation
}

# Function to create known latent states z
known.state.ms <- function(ms, notseen){
  # notseen: label for not seen
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(state)
}


# Function to reformat to m-array
marray <- function(ch, unobs = 0){
  ns <- length(table(ch)) - 1 + unobs
  no <- ncol(ch)
  out <- matrix(0, ncol = ns*(no-1)+1, nrow = ns*(no-1))
  # Remove capture histories of individuals that are marked at last occasion
  get.first <- function(x) min(which(x!=0))
  first <- apply(ch, 1, get.first)
  last.only <- which(first==no)
  if (length(last.only) > 0) ch <- ch[-last.only,]
  # Compute m-array
  for (i in 1:nrow(ch)){
    cap.occ <- which(ch[i,]!=0)
    state <- ch[i,cap.occ]
    if (length(state) == 1) {
      out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] <- out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] + 1
    }
    if (length(state) > 1) {
      for (t in 2:length(cap.occ)){
        out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*ns+state[t]] <- out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*ns+state[t]] + 1
      } # t
      if (max(cap.occ) < no){
        out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] <- out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] + 1
      } # if
    } # if
  } # t
  return(out)
}    

get.first <- function(x) min(which(x!=0))


#Creating a reordered matrix of distances between islands
distmat <- matrix(NA, nrow=10,ncol=10)
dist2 <- dist
dist2$SiteNo. <- c(10,7,9,8,6,3,5,4,1,2) #The numbers used for the islands in the model structure
dist2[11,]  <- c(NA,10,7,9,8,6,3,5,4,1,2,NA) #The numbers used for the islands in the model structure
for (i in 1:10){
  for (j in 1:10){
    inde <- which(as.integer(dist2[11,])==j)
    distmat[i,j] <- dist2[which(dist2$SiteNo.==i), inde]
  }}
rm("dist2")
########################
#######################








########################
########################
#Model 7
########################
#Estimating survival for 2-4-yr-olds jointly

#Creating an index matrix, to use in script
ind <- matrix(data=rep(1:10,10), ncol=10, nrow=10, byrow=T)
for (k in 1:10){ind[k,k] <- 11}

#Study characteristics
n.occasions <- 10 #Using the data from the last ten years only 
#(thus some birds will appear as they were first marked as 2/3/4-years-olds etc)
n.states <- 89 
n.obs <- 89 

#Removing the first years of data
enc0817 <- enc[, c(1:5, 19:31)] 

#Combining the numbers for adults as chicks as information about ringing stage is incorporated in the state
enc0817$No. <- enc0817$Ch + enc0817$Ad
#Selecting the relevant columns and changing the order
enc0817. <- enc0817[, c(18, 2:15, 19)] 

#Creating one encounter history for each individual
mat.enc <- matrix(NA, ncol = n.occasions, nrow = sum(enc0817.[,ncol(enc0817.)]))
i <- 1
for (l in 1:enc0817.[i,ncol(enc0817.)]){
  for (k in 1:n.occasions){
    mat.enc[l,k] <- enc0817.[i, (ncol(enc0817.)-n.occasions+k-1)]  } }

for (i in 2:nrow(enc0817.)){
  j <- 1 + sum(enc0817.[1:(i-1),ncol(enc0817.)])
  for (l in j:(sum(enc0817.[1:(i-1),ncol(enc0817.)])+enc0817.[i,ncol(enc0817.)])){
    for (k in 1:n.occasions){
      mat.enc[l,k] <- enc0817.[i, (ncol(enc0817.)-n.occasions+k-1)]    } }}

#Removing the rows with no events (only 0)
mat.enc2 <- mat.enc[which(rowSums(mat.enc)>0),]

# Compute vector with occasion of first capture
f <- apply(mat.enc2, 1, get.first)

#Checking the encounter histories
table(mat.enc2) #Need to check if there are any events that didn't happen, as some probabilities are low
length(table(mat.enc2))
#Reformating the data using the m-array function. 
ms.arr <- marray(mat.enc2, unobs=24) #Unobs needs to be no. of unobservable states + states not observed (despite observable)


p.all.years <- c(1,3,8) #Sites 1,3,8 should have p estimated for all years of the study

# Analysis of the model
# Specify model in BUGS language
sink("model7.jags") 
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi[fro, stage, t]: survival probability at age/state stage, having been at site fro, in period t
    # psi[fro, till]: movement probability from site fro to site till, only 10 of these as stage at next time-step is known and thereby constraints options to the 10 sites
    # p[...]: encounter probability. 
    # g: proportion of birds encountered (observed, captured or observed & captured) that where captured
    # -------------------------------------------------
    # States (S):
    # 1 alive at a1 s1
    # 2 alive at a2 s1
    # 3 alive at a3 s1
    # 4 alive at a4 s1
    # 5 alive at Ad_recent s1
    # 6 alive at Ad_previous s1
    # 7 alive at a1 s2
    # ...
    # 48 alive at Ad_previous s8
    # 49 alive & captured at a1 s9
    # 50 alive & captured at a2 s9
    # 51 alive & observed at a2 s9
    # 52 alive & observed & captured at a2 s9
    # 53 alive but not seen at a2
    # 54 alive & captured at a3 s10
    # ...
    # 85 alive Elsewhere at a2
    # 86 alive Elsewhere at a3
    # 87 alive Elsewhere at a4
    # 88 alive Elsewhere at Ad_previous
    # 89 dead
    #
    # Observations (O):
    # 1 seen at a1 s1
    # 2 seen at a2 s1
    # 3 seen at a3 s1
    # 4 seen at a4 s1
    # 5 seen at Ad s1
    # 6 seen at a1 s2
    # ...
    # 89/0 not seen
    # Unobservable states: 53,57,61,66, 71,75,79,84
    
    # -------------------------------------------------
    
    # Priors and constraints
    # Survival: random effects of time and island, additive effect of age
    for (fro in 1:10){
    for (stage in 1:6){ 
    for (t in 1:(n.occasions-1)){
    logit(phi[fro, stage, t]) <- mu.phi + eps.phi[fro, t] + del[stage] 
    }}}
    for (fro in 1:10){
    for (t in 1:(n.occasions-1)){
    eps.phi[fro, t] ~ dnorm(0, tau.phi)}}
    mean.phi ~ dunif(0, 1)                    # Prior for mean survival
    mu.phi <- log(mean.phi / (1-mean.phi))    # Logit transformation
    sigma.phi ~ dunif(0, 5)                       # Prior for standard deviation
    tau.phi <- pow(sigma.phi, -2)
    sigma2.phi <- pow(sigma.phi, 2)
    sigma2.phi.real <- sigma2.phi * pow(mean.phi, 2) * pow((1-mean.phi), 2) # Temporal variance on real scale
    
    del[1] <- 0                             # Priors for age/stage-effect
    for (stag in 2:4){ 
    del[stag] <- a.del[1]}
    a.del[1] ~ dnorm(0, 0.01)
    for (stag in 5:6){ 
    del[stag] <- a.del[stag-3]
    a.del[stag-3] ~ dnorm(0, 0.01)}
    for (fron in 1:10){  #Back-transform to the probability scale
    for (stage in 1:6){ 
    for (t in 1:(n.occasions-1)){
    phi.real[fron, stage, t] <- 1/(1+exp(-(mu.phi + eps.phi[fron,t] + del[stage])))    }}}
    for (st in 1:4){
    phiEl[st] ~ dunif(0, 1)}
    
    # Encounters at site 1-8: random effect of time and island, additive effect of age
    #Sites 1,3,8: 
    for (isl in p.all.years){ 
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[isl, stage, t]) <- mu.p + eps.p[isl, t] + delp[stage]     }}}
    for (isl in p.all.years){
    for (t in 1:(n.occasions-1)){
    eps.p[isl, t] ~ dnorm(0, tau.p)}}
    for (isln in p.all.years){  #Back-transform to the probability scale
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    p.real[isln, stage, t] <- 1/(1+exp(-(mu.p + eps.p[isln,t] + delp[stage])))    }}}
    #Site 2 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[2, stage, t]) <- mu.p + eps.p[2, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[2, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[2, stage, t] <- 1/(1+exp(-(mu.p + eps.p[2,t] + delp[stage])))    }}
    #Site 4 (should have p estimated for 2009-2011,2013-2017):
    for (stage in c(2:4,6)){ 
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    logit(p[4, stage, t]) <- mu.p + eps.p[4, t] + delp[stage]   }}
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    eps.p[4, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    p.real[4, stage, t] <- 1/(1+exp(-(mu.p + eps.p[4,t] + delp[stage])))    }}
    for (stage in c(2:4,6)){ 
    for (t in (n.occasions-6)){
    p[4, stage, t] <- 0   }}
    #Site 5 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[5, stage, t]) <- mu.p + eps.p[5, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[5, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[5, stage, t] <- 1/(1+exp(-(mu.p + eps.p[5,t] + delp[stage])))    }}
    #Site 6 (should have p estimated for 2011-2017):
    for (stage in c(2:4,6)){ 
    for (t in (n.occasions-7):(n.occasions-1)){
    logit(p[6, stage, t]) <- mu.p + eps.p[6, t] + delp[stage]   }}
    for (t in (n.occasions-7):(n.occasions-1)){
    eps.p[6, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in (n.occasions-7):(n.occasions-1)){
    p.real[6, stage, t] <- 1/(1+exp(-(mu.p + eps.p[6,t] + delp[stage])))    }}
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-8)){
    p[6, stage, t] <- 0   }}
    #Site 7 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[7, stage, t]) <- mu.p + eps.p[7, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[7, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[7, stage, t] <- 1/(1+exp(-(mu.p + eps.p[7,t] + delp[stage])))    }}
    
    mean.p ~ dunif(0, 1)                    # Prior for mean re-encounter
    mu.p <- log(mean.p / (1-mean.p))        # Logit transformation
    sigma.p ~ dunif(0, 5)                   # Prior for standard deviation
    tau.p <- pow(sigma.p, -2)
    sigma2.p <- pow(sigma.p, 2)
    sigma2.p.real <- sigma2.p * pow(mean.p, 2) * pow((1-mean.p), 2) # Temporal variance on real scale
    #delp[1] <- 0 #Re-encounters do not happen for age 1 [Shouldn't be necessary to define when using stage c(2:4,6) everywhere for p]
    #delp[5] <- 0 #Re-encounters do not happen for full-grown, just marked [Shouldn't be necessary to define when using stage c(2:4,6) everywhere for p]
    delp[2] <- 0 #Baseline age-class, which mu.p applies to
    for (stag in c(3,4)){ # Priors for age/stage-effect
    delp[stag] <- A.delp[stag-2]
    A.delp[stag-2] ~ dnorm(0, 0.01)}
    delp[6] <- A.delp[3]
    A.delp[3] ~ dnorm(0, 0.01)
    
    # Encounters at site 9-10:
    #Site 9 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p9[stage, t]) <- mu.p + eps.p9[t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p9[t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p9.real[stage, t] <- 1/(1+exp(-(mu.p + eps.p9[t] + delp[stage])))    }}
    #Site 10 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p10[stage, t]) <- mu.p + eps.p10[t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p10[t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p10.real[stage, t] <- 1/(1+exp(-(mu.p + eps.p10[t] + delp[stage])))    }}
    
    for (fd in 1:(n.occasions-5)){
    for (st in c(2:4,6)){
    g9[fd, st] <- 1  }}
    for (fd in (n.occasions-4):(n.occasions-1)){
    logit(g9[fd, 2]) <- t.g9[fd-n.occasions+5]
    logit(g9[fd, 3]) <- t.g9[fd-n.occasions+5] + a.g.3
    logit(g9[fd, 4]) <- t.g9[fd-n.occasions+5] + a.g.4
    logit(g9[fd, 6]) <- t.g9[fd-n.occasions+5] + a.g.6
    t.g9[fd-n.occasions+5] ~ dnorm(0, 0.01)
    g9.real[fd,2] <- 1/(1+exp(-(t.g9[fd-n.occasions+5]))) 
    g9.real[fd,3] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.3)))
    g9.real[fd,4] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.4)))
    g9.real[fd,6] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.6)))}
    for (fd in 1:(n.occasions-4)){ 
    for (st in c(2:4,6)){
    g10[fd, st] <- 1  }}
    logit(g10[n.occasions-3, 2]) <- t.g10
    logit(g10[n.occasions-3, 3]) <- t.g10 + a.g.3
    logit(g10[n.occasions-3, 4]) <- t.g10 + a.g.4
    logit(g10[n.occasions-3, 6]) <- t.g10 + a.g.6
    t.g10 ~ dnorm(0, 0.01)
    g10.real[n.occasions-3,2] <- 1/(1+exp(-(t.g10))) 
    g10.real[n.occasions-3,3] <- 1/(1+exp(-(t.g10 + a.g.3)))
    g10.real[n.occasions-3,4] <- 1/(1+exp(-(t.g10 + a.g.4)))
    g10.real[n.occasions-3,6] <- 1/(1+exp(-(t.g10 + a.g.6)))
    for (fd in (n.occasions-2):(n.occasions-1)){ 
    for (st in c(2:4,6)){
    g10[fd, st] <- 1  }}
    a.g.3 ~ dnorm(0, 0.01)
    a.g.4 ~ dnorm(0, 0.01)
    a.g.6 ~ dnorm(0, 0.01)
    
    # Encounters at site 11:
    for (st in 1:4){ 
    pEl[st] ~ dunif(0, 1) }
    
    # Transitions: multinomial logit
    # Normal priors on logit of all but one transition probs
    for (fro in 1:10){
    for (till in 1:(fro-1)){
    lpsi1[fro, till] <- repuls[fro] + attract[till]
    lpsi2[fro, till] <- theta2 + repuls[fro] + attract[till]
    lpsi3[fro, till] <- theta3 + repuls[fro] + attract[till]
    lpsi4[fro, till] <- theta4 + repuls[fro] + attract[till]
    lpsiAa[fro, till] <- thetaAa + repuls[fro] + attract[till]
    lpsiAb[fro, till] <- thetaAb + repuls[fro] + attract[till]}
    for (till in (fro+1):10){
    lpsi1[fro, till] <- repuls[fro] + attract[till]
    lpsi2[fro, till] <- theta2 + repuls[fro] + attract[till]
    lpsi3[fro, till] <- theta3 + repuls[fro] + attract[till]
    lpsi4[fro, till] <- theta4 + repuls[fro] + attract[till]
    lpsiAa[fro, till] <- thetaAa + repuls[fro] + attract[till]
    lpsiAb[fro, till] <- thetaAb + repuls[fro] + attract[till]} 
    lpsi1[fro, 11] <- repuls[fro] + attract[11]
    lpsi2[fro, 11] <- theta2 + repuls[fro] + attract[11]
    lpsi3[fro, 11] <- theta3 + repuls[fro] + attract[11]
    lpsi4[fro, 11] <- theta4 + repuls[fro] + attract[11]
    lpsiAa[fro, 11] <- thetaAa + repuls[fro] + attract[11]
    lpsiAb[fro, 11] <- thetaAb + repuls[fro] + attract[11]}
    
    for (fro in 1:10){
    for (till in 1:(fro-1)){
    psi1[fro, till] <- exp(lpsi1[fro, till]) / (1 + exp(lpsi1[fro,ind[fro,1]]) + exp(lpsi1[fro,ind[fro,2]]) + exp(lpsi1[fro,ind[fro,3]]) + exp(lpsi1[fro,ind[fro,4]]) + exp(lpsi1[fro,ind[fro,5]]) + exp(lpsi1[fro,ind[fro,6]]) + exp(lpsi1[fro,ind[fro,7]]) + exp(lpsi1[fro,ind[fro,8]]) + exp(lpsi1[fro,ind[fro,9]]) + exp(lpsi1[fro,ind[fro,10]]))
    psi2[fro, till] <- exp(lpsi2[fro, till]) / (1 + exp(lpsi2[fro,ind[fro,1]]) + exp(lpsi2[fro,ind[fro,2]]) + exp(lpsi2[fro,ind[fro,3]]) + exp(lpsi2[fro,ind[fro,4]]) + exp(lpsi2[fro,ind[fro,5]]) + exp(lpsi2[fro,ind[fro,6]]) + exp(lpsi2[fro,ind[fro,7]]) + exp(lpsi2[fro,ind[fro,8]]) + exp(lpsi2[fro,ind[fro,9]]) + exp(lpsi2[fro,ind[fro,10]]))
    psi3[fro, till] <- exp(lpsi3[fro, till]) / (1 + exp(lpsi3[fro,ind[fro,1]]) + exp(lpsi3[fro,ind[fro,2]]) + exp(lpsi3[fro,ind[fro,3]]) + exp(lpsi3[fro,ind[fro,4]]) + exp(lpsi3[fro,ind[fro,5]]) + exp(lpsi3[fro,ind[fro,6]]) + exp(lpsi3[fro,ind[fro,7]]) + exp(lpsi3[fro,ind[fro,8]]) + exp(lpsi3[fro,ind[fro,9]]) + exp(lpsi3[fro,ind[fro,10]]))
    psi4[fro, till] <- exp(lpsi4[fro, till]) / (1 + exp(lpsi4[fro,ind[fro,1]]) + exp(lpsi4[fro,ind[fro,2]]) + exp(lpsi4[fro,ind[fro,3]]) + exp(lpsi4[fro,ind[fro,4]]) + exp(lpsi4[fro,ind[fro,5]]) + exp(lpsi4[fro,ind[fro,6]]) + exp(lpsi4[fro,ind[fro,7]]) + exp(lpsi4[fro,ind[fro,8]]) + exp(lpsi4[fro,ind[fro,9]]) + exp(lpsi4[fro,ind[fro,10]]))
    psiAa[fro, till] <- exp(lpsiAa[fro, till]) / (1 + exp(lpsiAa[fro,ind[fro,1]]) + exp(lpsiAa[fro,ind[fro,2]]) + exp(lpsiAa[fro,ind[fro,3]]) + exp(lpsiAa[fro,ind[fro,4]]) + exp(lpsiAa[fro,ind[fro,5]]) + exp(lpsiAa[fro,ind[fro,6]]) + exp(lpsiAa[fro,ind[fro,7]]) + exp(lpsiAa[fro,ind[fro,8]]) + exp(lpsiAa[fro,ind[fro,9]]) + exp(lpsiAa[fro,ind[fro,10]]))
    psiAb[fro, till] <- exp(lpsiAb[fro, till]) / (1 + exp(lpsiAb[fro,ind[fro,1]]) + exp(lpsiAb[fro,ind[fro,2]]) + exp(lpsiAb[fro,ind[fro,3]]) + exp(lpsiAb[fro,ind[fro,4]]) + exp(lpsiAb[fro,ind[fro,5]]) + exp(lpsiAb[fro,ind[fro,6]]) + exp(lpsiAb[fro,ind[fro,7]]) + exp(lpsiAb[fro,ind[fro,8]]) + exp(lpsiAb[fro,ind[fro,9]]) + exp(lpsiAb[fro,ind[fro,10]]))  }
    for (till in (fro+1):11){
    psi1[fro, till] <- exp(lpsi1[fro, till]) / (1 + exp(lpsi1[fro,ind[fro,1]]) + exp(lpsi1[fro,ind[fro,2]]) + exp(lpsi1[fro,ind[fro,3]]) + exp(lpsi1[fro,ind[fro,4]]) + exp(lpsi1[fro,ind[fro,5]]) + exp(lpsi1[fro,ind[fro,6]]) + exp(lpsi1[fro,ind[fro,7]]) + exp(lpsi1[fro,ind[fro,8]]) + exp(lpsi1[fro,ind[fro,9]]) + exp(lpsi1[fro,ind[fro,10]]))
    psi2[fro, till] <- exp(lpsi2[fro, till]) / (1 + exp(lpsi2[fro,ind[fro,1]]) + exp(lpsi2[fro,ind[fro,2]]) + exp(lpsi2[fro,ind[fro,3]]) + exp(lpsi2[fro,ind[fro,4]]) + exp(lpsi2[fro,ind[fro,5]]) + exp(lpsi2[fro,ind[fro,6]]) + exp(lpsi2[fro,ind[fro,7]]) + exp(lpsi2[fro,ind[fro,8]]) + exp(lpsi2[fro,ind[fro,9]]) + exp(lpsi2[fro,ind[fro,10]]))
    psi3[fro, till] <- exp(lpsi3[fro, till]) / (1 + exp(lpsi3[fro,ind[fro,1]]) + exp(lpsi3[fro,ind[fro,2]]) + exp(lpsi3[fro,ind[fro,3]]) + exp(lpsi3[fro,ind[fro,4]]) + exp(lpsi3[fro,ind[fro,5]]) + exp(lpsi3[fro,ind[fro,6]]) + exp(lpsi3[fro,ind[fro,7]]) + exp(lpsi3[fro,ind[fro,8]]) + exp(lpsi3[fro,ind[fro,9]]) + exp(lpsi3[fro,ind[fro,10]]))
    psi4[fro, till] <- exp(lpsi4[fro, till]) / (1 + exp(lpsi4[fro,ind[fro,1]]) + exp(lpsi4[fro,ind[fro,2]]) + exp(lpsi4[fro,ind[fro,3]]) + exp(lpsi4[fro,ind[fro,4]]) + exp(lpsi4[fro,ind[fro,5]]) + exp(lpsi4[fro,ind[fro,6]]) + exp(lpsi4[fro,ind[fro,7]]) + exp(lpsi4[fro,ind[fro,8]]) + exp(lpsi4[fro,ind[fro,9]]) + exp(lpsi4[fro,ind[fro,10]]))
    psiAa[fro, till] <- exp(lpsiAa[fro, till]) / (1 + exp(lpsiAa[fro,ind[fro,1]]) + exp(lpsiAa[fro,ind[fro,2]]) + exp(lpsiAa[fro,ind[fro,3]]) + exp(lpsiAa[fro,ind[fro,4]]) + exp(lpsiAa[fro,ind[fro,5]]) + exp(lpsiAa[fro,ind[fro,6]]) + exp(lpsiAa[fro,ind[fro,7]]) + exp(lpsiAa[fro,ind[fro,8]]) + exp(lpsiAa[fro,ind[fro,9]]) + exp(lpsiAa[fro,ind[fro,10]]))
    psiAb[fro, till] <- exp(lpsiAb[fro, till]) / (1 + exp(lpsiAb[fro,ind[fro,1]]) + exp(lpsiAb[fro,ind[fro,2]]) + exp(lpsiAb[fro,ind[fro,3]]) + exp(lpsiAb[fro,ind[fro,4]]) + exp(lpsiAb[fro,ind[fro,5]]) + exp(lpsiAb[fro,ind[fro,6]]) + exp(lpsiAb[fro,ind[fro,7]]) + exp(lpsiAb[fro,ind[fro,8]]) + exp(lpsiAb[fro,ind[fro,9]]) + exp(lpsiAb[fro,ind[fro,10]]))  }
    #The case below is when fro=till, but rewritten to match JAGS syntax
    psi1[fro,fro] <- 1- sum(psi1[fro,ind[fro,]]) 
    psi2[fro,fro] <- 1- sum(psi2[fro,ind[fro,]]) 
    psi3[fro,fro] <- 1- sum(psi3[fro,ind[fro,]]) 
    psi4[fro,fro] <- 1- sum(psi4[fro,ind[fro,]]) 
    psiAa[fro,fro] <- 1- sum(psiAa[fro,ind[fro,]])   
    psiAb[fro,fro] <- 1- sum(psiAb[fro,ind[fro,]])   
    } 
    
    #Priors for transitions
    theta2 ~ dnorm(0, 0.01)               # Prior for age-effect in psi
    theta3 ~ dnorm(0, 0.01)  
    theta4 ~ dnorm(0, 0.01)  
    thetaAa ~ dnorm(0, 0.01)  
    thetaAb ~ dnorm(0, 0.01) 
    #mean.psi ~ dunif(0, 1)                # Prior for mean movement prob - not used as no transition makes sense to set as baseline
    #mu <- log(mean.psi / (1-mean.psi))     # Logit transformation
    for (w in 1:10){
    repuls[w] ~ dnorm(0, 0.01) 
    attract[w] ~ dnorm(0, 0.01) 
    }
    attract[11] ~ dnorm(0, 0.01)
    
    #Transitions from Elsewhere
    for (till in 1:10){
    lpsi2El[till] <- muEl 
    lpsi3El[till] <- muEl + theta3
    lpsi4El[till] <- muEl + theta4
    lpsiAbEl[till] <- muEl + thetaAb
    }
    for (till in 1:10){
    psi2El[till] <- exp(lpsi2El[till]) / (1 + exp(lpsi2El[1]) + exp(lpsi2El[2]) + exp(lpsi2El[3]) + exp(lpsi2El[4]) + exp(lpsi2El[5]) + exp(lpsi2El[6]) + exp(lpsi2El[7]) + exp(lpsi2El[8]) + exp(lpsi2El[9]) + exp(lpsi2El[10]))
    psi3El[till] <- exp(lpsi3El[till]) / (1 + exp(lpsi3El[1]) + exp(lpsi3El[2]) + exp(lpsi3El[3]) + exp(lpsi3El[4]) + exp(lpsi3El[5]) + exp(lpsi3El[6]) + exp(lpsi3El[7]) + exp(lpsi3El[8]) + exp(lpsi3El[9]) + exp(lpsi3El[10]))
    psi4El[till] <- exp(lpsi4El[till]) / (1 + exp(lpsi4El[1]) + exp(lpsi4El[2]) + exp(lpsi4El[3]) + exp(lpsi4El[4]) + exp(lpsi4El[5]) + exp(lpsi4El[6]) + exp(lpsi4El[7]) + exp(lpsi4El[8]) + exp(lpsi4El[9]) + exp(lpsi4El[10]))
    psiAbEl[till] <- exp(lpsiAbEl[till]) / (1 + exp(lpsiAbEl[1]) + exp(lpsiAbEl[2]) + exp(lpsiAbEl[3]) + exp(lpsiAbEl[4]) + exp(lpsiAbEl[5]) + exp(lpsiAbEl[6]) + exp(lpsiAbEl[7]) + exp(lpsiAbEl[8]) + exp(lpsiAbEl[9]) + exp(lpsiAbEl[10]))
    }
    # Calculate the last transition probability
    psi2El[11] <- 1- sum(psi2El[1:10]) 
    psi3El[11] <- 1- sum(psi3El[1:10]) 
    psi4El[11] <- 1- sum(psi4El[1:10]) 
    psiAbEl[11] <- 1- sum(psiAbEl[1:10])   
    #Prior for transitions from Elsewhere    
    mean.psiEl ~ dunif(0, 1)                   # Prior for mean movement prob from Elsewhere
    muEl <- log(mean.psiEl / (1-mean.psiEl))     # Logit transformation
    
    
    # Define state-transition and observation matrices 	
    # Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
    for (d in 1:88){
    for (m in 1:8){
    psi[d,t,1+(6*(m-1))] <- 0
    psi[d,t,5+(6*(m-1))] <- 0
    } #m
    psi[d,t,49] <- 0
    psi[d,t,62] <- 0
    psi[d,t,67] <- 0
    psi[d,t,80] <- 0
    } #d
    
    for (sf in 1:8){  
    for (s in 1:8){
    psi[(6*(sf-1))+1,t,(s*6-4)] <- phi[sf,1,t]*psi1[sf,s]
    psi[(6*(sf-1))+1,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+1,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+1,t,(s*6)] <- 0
    psi[(6*(sf-1))+2,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+2,t,(s*6-3)] <- phi[sf,2,t]*psi2[sf,s]
    psi[(6*(sf-1))+2,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+2,t,(s*6)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-2)] <- phi[sf,3,t]*psi3[sf,s]
    psi[(6*(sf-1))+3,t,(s*6)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+4,t,(s*6)] <- phi[sf,4,t]*psi4[sf,s]
    psi[(6*(sf-1))+5,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+5,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+5,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+5,t,(s*6)] <- phi[sf,5,t]*psiAa[sf,s]
    psi[(6*(sf-1))+6,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+6,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+6,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+6,t,(s*6)] <- phi[sf,6,t]*psiAb[sf,s]
    } #s
    
    psi[(6*(sf-1))+1,t,50] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*g9[t, 2]
    psi[(6*(sf-1))+1,t,51] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[(6*(sf-1))+1,t,52] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[(6*(sf-1))+1,t,53] <- phi[sf,1,t]*psi1[sf,9]*(1-p9[2,t])
    for (ex in c(54:61,63:66)){
    psi[(6*(sf-1))+1,t,ex] <- 0  }
    psi[(6*(sf-1))+1,t,68] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*g10[t, 2]
    psi[(6*(sf-1))+1,t,69] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[(6*(sf-1))+1,t,70] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[(6*(sf-1))+1,t,71] <- phi[sf,1,t]*psi1[sf,10]*(1-p10[2,t])
    for (ex in c(72:79,81:84)){
    psi[(6*(sf-1))+1,t,ex] <- 0  }
    psi[(6*(sf-1))+1,t,85] <- phi[sf,1,t]*psi1[sf,11]
    psi[(6*(sf-1))+1,t,86] <- 0
    psi[(6*(sf-1))+1,t,87] <- 0
    psi[(6*(sf-1))+1,t,88] <- 0
    
    for (ew in 50:53){
    psi[(6*(sf-1))+2,t,ew] <- 0}
    psi[(6*(sf-1))+2,t,54] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*g9[t, 3]
    psi[(6*(sf-1))+2,t,55] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[(6*(sf-1))+2,t,56] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[(6*(sf-1))+2,t,57] <- phi[sf,2,t]*psi2[sf,9]*(1-p9[3,t])
    for (er in c(58:61,63:66)){
    psi[(6*(sf-1))+2,t,er] <- 0  }
    for (ew in 68:71){
    psi[(6*(sf-1))+2,t,ew] <- 0}
    psi[(6*(sf-1))+2,t,72] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*g10[t, 3]
    psi[(6*(sf-1))+2,t,73] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[(6*(sf-1))+2,t,74] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[(6*(sf-1))+2,t,75] <- phi[sf,2,t]*psi2[sf,10]*(1-p10[3,t])
    for (er in c(76:79,81:84)){
    psi[(6*(sf-1))+2,t,er] <- 0  }
    psi[(6*(sf-1))+2,t,85] <- 0
    psi[(6*(sf-1))+2,t,86] <- phi[sf,2,t]*psi2[sf,11]
    psi[(6*(sf-1))+2,t,87] <- 0
    psi[(6*(sf-1))+2,t,88] <- 0
    
    for (ed in 50:57){
    psi[(6*(sf-1))+3,t,ed] <- 0}
    psi[(6*(sf-1))+3,t,58] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*g9[t, 4]
    psi[(6*(sf-1))+3,t,59] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[(6*(sf-1))+3,t,60] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[(6*(sf-1))+3,t,61] <- phi[sf,3,t]*psi3[sf,9]*(1-p9[4,t])
    for (ef in 63:66){
    psi[(6*(sf-1))+3,t,ef] <- 0 }
    for (ed in 68:75){
    psi[(6*(sf-1))+3,t,ed] <- 0}
    psi[(6*(sf-1))+3,t,76] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*g10[t, 4]
    psi[(6*(sf-1))+3,t,77] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[(6*(sf-1))+3,t,78] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[(6*(sf-1))+3,t,79] <- phi[sf,3,t]*psi3[sf,10]*(1-p10[4,t])
    for (ef in 81:84){
    psi[(6*(sf-1))+3,t,ef] <- 0 }
    psi[(6*(sf-1))+3,t,85] <- 0
    psi[(6*(sf-1))+3,t,86] <- 0
    psi[(6*(sf-1))+3,t,87] <- phi[sf,3,t]*psi3[sf,11]
    psi[(6*(sf-1))+3,t,88] <- 0
    
    for (es in 50:61){
    psi[(6*(sf-1))+4,t,es] <- 0}
    psi[(6*(sf-1))+4,t,63] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+4,t,64] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+4,t,65] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+4,t,66] <- phi[sf,4,t]*psi4[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+4,t,es] <- 0}
    psi[(6*(sf-1))+4,t,81] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+4,t,82] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+4,t,83] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+4,t,84] <- phi[sf,4,t]*psi4[sf,10]*(1-p10[6,t])
    psi[(6*(sf-1))+4,t,85] <- 0
    psi[(6*(sf-1))+4,t,86] <- 0
    psi[(6*(sf-1))+4,t,87] <- 0
    psi[(6*(sf-1))+4,t,88] <- phi[sf,4,t]*psi4[sf,11]
    
    for (es in 50:61){
    psi[(6*(sf-1))+5,t,es] <- 0}
    psi[(6*(sf-1))+5,t,63] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+5,t,64] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+5,t,65] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+5,t,66] <- phi[sf,5,t]*psiAa[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+5,t,es] <- 0}
    psi[(6*(sf-1))+5,t,81] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+5,t,82] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+5,t,83] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+5,t,84] <- phi[sf,5,t]*psiAa[sf,10]*(1-p10[6,t])
    psi[(6*(sf-1))+5,t,85] <- 0
    psi[(6*(sf-1))+5,t,86] <- 0
    psi[(6*(sf-1))+5,t,87] <- 0
    psi[(6*(sf-1))+5,t,88] <- phi[sf,5,t]*psiAa[sf,11]
    
    for (es in 50:61){
    psi[(6*(sf-1))+6,t,es] <- 0}
    psi[(6*(sf-1))+6,t,63] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+6,t,64] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+6,t,65] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+6,t,66] <- phi[sf,6,t]*psiAb[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+6,t,es] <- 0}
    psi[(6*(sf-1))+6,t,81] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+6,t,82] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+6,t,83] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+6,t,84] <- phi[sf,6,t]*psiAb[sf,10]*(1-p10[6,t])
    psi[(6*(sf-1))+6,t,85] <- 0
    psi[(6*(sf-1))+6,t,86] <- 0
    psi[(6*(sf-1))+6,t,87] <- 0
    psi[(6*(sf-1))+6,t,88] <- phi[sf,6,t]*psiAb[sf,11]
    } #sf
    
    
    for (s in 1:8){
    psi[49,t,(s*6-4)] <- phi[9,1,t]*psi1[9,s]
    psi[49,t,(s*6-3)] <- 0
    psi[49,t,(s*6-2)] <- 0
    psi[49,t,(s*6)] <- 0
    for (i in 1:4){
    psi[49+i,t,(s*6-4)] <- 0
    psi[49+i,t,(s*6-3)] <- phi[9,2,t]*psi2[9,s]
    psi[49+i,t,(s*6-2)] <- 0
    psi[49+i,t,(s*6)] <- 0
    psi[53+i,t,(s*6-4)] <- 0
    psi[53+i,t,(s*6-3)] <- 0
    psi[53+i,t,(s*6-2)] <- phi[9,3,t]*psi3[9,s]
    psi[53+i,t,(s*6)] <- 0
    psi[57+i,t,(s*6-4)] <- 0
    psi[57+i,t,(s*6-3)] <- 0
    psi[57+i,t,(s*6-2)] <- 0
    psi[57+i,t,(s*6)] <- phi[9,4,t]*psi4[9,s] 
    psi[62+i,t,(s*6-4)] <- 0
    psi[62+i,t,(s*6-3)] <- 0
    psi[62+i,t,(s*6-2)] <- 0
    psi[62+i,t,(s*6)] <- phi[9,6,t]*psiAb[9,s]
    } #i
    psi[62,t,(s*6-4)] <- 0
    psi[62,t,(s*6-3)] <- 0
    psi[62,t,(s*6-2)] <- 0
    psi[62,t,(s*6)] <- phi[9,5,t]*psiAa[9,s]
    } #s
    
    for (s in 1:8){
    psi[67,t,(s*6-4)] <- phi[10,1,t]*psi1[10,s]
    psi[67,t,(s*6-3)] <- 0
    psi[67,t,(s*6-2)] <- 0
    psi[67,t,(s*6)] <- 0
    for (i in 1:4){
    psi[67+i,t,(s*6-4)] <- 0
    psi[67+i,t,(s*6-3)] <- phi[10,2,t]*psi2[10,s]
    psi[67+i,t,(s*6-2)] <- 0
    psi[67+i,t,(s*6)] <- 0
    psi[71+i,t,(s*6-4)] <- 0
    psi[71+i,t,(s*6-3)] <- 0
    psi[71+i,t,(s*6-2)] <- phi[10,3,t]*psi3[10,s]
    psi[71+i,t,(s*6)] <- 0
    psi[75+i,t,(s*6-4)] <- 0
    psi[75+i,t,(s*6-3)] <- 0
    psi[75+i,t,(s*6-2)] <- 0
    psi[75+i,t,(s*6)] <- phi[10,4,t]*psi4[10,s] 
    psi[80+i,t,(s*6-4)] <- 0
    psi[80+i,t,(s*6-3)] <- 0
    psi[80+i,t,(s*6-2)] <- 0
    psi[80+i,t,(s*6)] <- phi[10,6,t]*psiAb[10,s]
    } #i
    psi[80,t,(s*6-4)] <- 0
    psi[80,t,(s*6-3)] <- 0
    psi[80,t,(s*6-2)] <- 0
    psi[80,t,(s*6)] <- phi[10,5,t]*psiAa[10,s]
    } #s
    
    psi[49,t,50] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*g9[t, 2]
    psi[49,t,51] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[49,t,52] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[49,t,53] <- phi[9,1,t]*psi1[9,9]*(1-p9[2,t])
    for (ex in c(54:61,63:66)){
    psi[49,t,ex] <- 0  }
    
    for (i in 1:4){
    for (ew in 50:53){
    psi[49+i,t,ew] <- 0  }
    psi[49+i,t,54] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*g9[t, 3]
    psi[49+i,t,55] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[49+i,t,56] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[49+i,t,57] <- phi[9,2,t]*psi2[9,9]*(1-p9[3,t])
    for (er in c(58:61,63:66)){
    psi[49+i,t,er] <- 0  }
    for (ed in 50:57){
    psi[53+i,t,ed] <- 0}
    psi[53+i,t,58] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*g9[t, 4]
    psi[53+i,t,59] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[53+i,t,60] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[53+i,t,61] <- phi[9,3,t]*psi3[9,9]*(1-p9[4,t])
    for (ef in 63:66){
    psi[53+i,t,ef] <- 0 }
    for (es in 50:61){
    psi[57+i,t,es] <- 0}
    psi[57+i,t,63] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*g9[t, 6]
    psi[57+i,t,64] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[57+i,t,65] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[57+i,t,66] <- phi[9,4,t]*psi4[9,9]*(1-p9[6,t])
    for (es in 50:61){
    psi[62+i,t,es] <- 0}
    psi[62+i,t,63] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*g9[t, 6]
    psi[62+i,t,64] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62+i,t,65] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62+i,t,66] <- phi[9,6,t]*psiAb[9,9]*(1-p9[6,t])
    } #i
    
    for (es in 50:61){
    psi[62,t,es] <- 0}
    psi[62,t,63] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*g9[t, 6]
    psi[62,t,64] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62,t,65] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62,t,66] <- phi[9,5,t]*psiAa[9,9]*(1-p9[6,t])
    
    
    psi[49,t,50+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*g10[t, 2]
    psi[49,t,51+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[49,t,52+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[49,t,53+18] <- phi[9,1,t]*psi1[9,10]*(1-p10[2,t])
    for (ex in c(54:61,63:66)){
    psi[49,t,ex+18] <- 0  }
    
    for (i in 1:4){
    for (ew in 50:53){
    psi[49+i,t,ew+18] <- 0  }
    psi[49+i,t,54+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*g10[t, 3]
    psi[49+i,t,55+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[49+i,t,56+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[49+i,t,57+18] <- phi[9,2,t]*psi2[9,10]*(1-p10[3,t])
    for (er in c(58:61,63:66)){
    psi[49+i,t,er+18] <- 0  }
    for (ed in 50:57){
    psi[53+i,t,ed+18] <- 0}
    psi[53+i,t,58+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*g10[t, 4]
    psi[53+i,t,59+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[53+i,t,60+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[53+i,t,61+18] <- phi[9,3,t]*psi3[9,10]*(1-p10[4,t])
    for (ef in 63:66){
    psi[53+i,t,ef+18] <- 0 }
    for (es in 50:61){
    psi[57+i,t,es+18] <- 0}
    psi[57+i,t,63+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*g10[t, 6]
    psi[57+i,t,64+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[57+i,t,65+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[57+i,t,66+18] <- phi[9,4,t]*psi4[9,10]*(1-p10[6,t])
    for (es in 50:61){
    psi[62+i,t,es+18] <- 0}
    psi[62+i,t,63+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*g10[t, 6]
    psi[62+i,t,64+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62+i,t,65+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62+i,t,66+18] <- phi[9,6,t]*psiAb[9,10]*(1-p10[6,t])
    } #i
    
    for (es in 50:61){
    psi[62,t,es+18] <- 0}
    psi[62,t,63+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*g10[t, 6]
    psi[62,t,64+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62,t,65+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62,t,66+18] <- phi[9,5,t]*psiAa[9,10]*(1-p10[6,t])
    
    
    psi[49,t,(89-4)] <- phi[9,1,t]*psi1[9,11]
    psi[49,t,(89-3)] <- 0
    psi[49,t,(89-2)] <- 0
    psi[49,t,(88)] <- 0
    for (i in 1:4){
    psi[49+i,t,(89-4)] <- 0
    psi[49+i,t,(89-3)] <- phi[9,2,t]*psi2[9,11]
    psi[49+i,t,(89-2)] <- 0
    psi[49+i,t,(88)] <- 0
    psi[53+i,t,(89-4)] <- 0
    psi[53+i,t,(89-3)] <- 0
    psi[53+i,t,(89-2)] <- phi[9,3,t]*psi3[9,11]
    psi[53+i,t,(88)] <- 0
    psi[57+i,t,(89-4)] <- 0
    psi[57+i,t,(89-3)] <- 0
    psi[57+i,t,(89-2)] <- 0
    psi[57+i,t,(88)] <- phi[9,4,t]*psi4[9,11] 
    psi[62+i,t,(89-4)] <- 0
    psi[62+i,t,(89-3)] <- 0
    psi[62+i,t,(89-2)] <- 0
    psi[62+i,t,(88)] <- phi[9,6,t]*psiAb[9,11]
    } #i
    psi[62,t,(89-4)] <- 0
    psi[62,t,(89-3)] <- 0
    psi[62,t,(89-2)] <- 0
    psi[62,t,(88)] <- phi[9,5,t]*psiAa[9,11]
    
    
    psi[67,t,50] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*g9[t, 2]
    psi[67,t,51] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[67,t,52] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[67,t,53] <- phi[10,1,t]*psi1[10,9]*(1-p9[2,t])
    for (ex in c(72:79,81:84)){
    psi[67,t,ex-18] <- 0  }
    
    for (i in 1:4){
    for (ew in 68:71){
    psi[67+i,t,ew-18] <- 0  }
    psi[67+i,t,72-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*g9[t, 3]
    psi[67+i,t,73-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[67+i,t,74-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[67+i,t,75-18] <- phi[10,2,t]*psi2[10,9]*(1-p9[3,t])
    for (er in c(76:79,81:84)){
    psi[67+i,t,er-18] <- 0  }
    for (ed in 68:75){
    psi[71+i,t,ed-18] <- 0}
    psi[71+i,t,76-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*g9[t, 4]
    psi[71+i,t,77-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[71+i,t,78-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[71+i,t,79-18] <- phi[10,3,t]*psi3[10,9]*(1-p9[4,t])
    for (ef in 81:84){
    psi[71+i,t,ef-18] <- 0 }
    for (es in 68:79){
    psi[75+i,t,es-18] <- 0}
    psi[75+i,t,81-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*g9[t, 6]
    psi[75+i,t,82-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[75+i,t,83-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[75+i,t,84-18] <- phi[10,4,t]*psi4[10,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[80+i,t,es-18] <- 0}
    psi[80+i,t,81-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*g9[t, 6]
    psi[80+i,t,82-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80+i,t,83-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80+i,t,84-18] <- phi[10,6,t]*psiAb[10,9]*(1-p9[6,t])
    } #i
    
    for (es in 68:79){
    psi[80,t,es-18] <- 0}
    psi[80,t,81-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*g9[t, 6]
    psi[80,t,82-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80,t,83-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80,t,84-18] <- phi[10,5,t]*psiAa[10,9]*(1-p9[6,t])
    
    
    psi[67,t,68] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*g10[t, 2]
    psi[67,t,69] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[67,t,70] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[67,t,71] <- phi[10,1,t]*psi1[10,10]*(1-p10[2,t])
    for (ex in c(72:79,81:84)){
    psi[67,t,ex] <- 0  }
    
    for (i in 1:4){
    for (ew in 68:71){
    psi[67+i,t,ew] <- 0  }
    psi[67+i,t,72] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*g10[t, 3]
    psi[67+i,t,73] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[67+i,t,74] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[67+i,t,75] <- phi[10,2,t]*psi2[10,10]*(1-p10[3,t])
    for (er in c(76:79,81:84)){
    psi[67+i,t,er] <- 0  }
    for (ed in 68:75){
    psi[71+i,t,ed] <- 0}
    psi[71+i,t,76] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*g10[t, 4]
    psi[71+i,t,77] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[71+i,t,78] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[71+i,t,79] <- phi[10,3,t]*psi3[10,10]*(1-p10[4,t])
    for (ef in 81:84){
    psi[71+i,t,ef] <- 0 }
    for (es in 68:79){
    psi[75+i,t,es] <- 0}
    psi[75+i,t,81] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*g10[t, 6]
    psi[75+i,t,82] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[75+i,t,83] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[75+i,t,84] <- phi[10,4,t]*psi4[10,10]*(1-p10[6,t])
    for (es in 68:79){
    psi[80+i,t,es] <- 0}
    psi[80+i,t,81] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*g10[t, 6]
    psi[80+i,t,82] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80+i,t,83] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80+i,t,84] <- phi[10,6,t]*psiAb[10,10]*(1-p10[6,t])
    } #i
    
    for (es in 68:79){
    psi[80,t,es] <- 0}
    psi[80,t,81] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*g10[t, 6]
    psi[80,t,82] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80,t,83] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80,t,84] <- phi[10,5,t]*psiAa[10,10]*(1-p10[6,t])
    
    
    psi[67,t,(89-4)] <- phi[10,1,t]*psi1[10,11]
    psi[67,t,(89-3)] <- 0
    psi[67,t,(89-2)] <- 0
    psi[67,t,(88)] <- 0
    for (i in 1:4){
    psi[67+i,t,(89-4)] <- 0
    psi[67+i,t,(89-3)] <- phi[10,2,t]*psi2[10,11]
    psi[67+i,t,(89-2)] <- 0
    psi[67+i,t,(88)] <- 0
    psi[71+i,t,(89-4)] <- 0
    psi[71+i,t,(89-3)] <- 0
    psi[71+i,t,(89-2)] <- phi[10,3,t]*psi3[10,11]
    psi[71+i,t,(88)] <- 0
    psi[75+i,t,(89-4)] <- 0
    psi[75+i,t,(89-3)] <- 0
    psi[75+i,t,(89-2)] <- 0
    psi[75+i,t,(88)] <- phi[10,4,t]*psi4[10,11] 
    psi[80+i,t,(89-4)] <- 0
    psi[80+i,t,(89-3)] <- 0
    psi[80+i,t,(89-2)] <- 0
    psi[80+i,t,(88)] <- phi[10,6,t]*psiAb[10,11]
    } #i
    psi[80,t,(89-4)] <- 0
    psi[80,t,(89-3)] <- 0
    psi[80,t,(89-2)] <- 0
    psi[80,t,(88)] <- phi[10,5,t]*psiAa[10,11]
    
    for (s in 1:8){
    psi[85,t,(s*6-4)] <- 0
    psi[85,t,(s*6-3)] <- phiEl[1]*psi2El[s]
    psi[85,t,(s*6-2)] <- 0
    psi[85,t,(s*6)] <- 0
    psi[86,t,(s*6-4)] <- 0
    psi[86,t,(s*6-3)] <- 0
    psi[86,t,(s*6-2)] <- phiEl[2]*psi3El[s]
    psi[86,t,(s*6)] <- 0
    psi[87,t,(s*6-4)] <- 0
    psi[87,t,(s*6-3)] <- 0
    psi[87,t,(s*6-2)] <- 0
    psi[87,t,(s*6)] <- phiEl[3]*psi4El[s]
    psi[88,t,(s*6-4)] <- 0
    psi[88,t,(s*6-3)] <- 0
    psi[88,t,(s*6-2)] <- 0
    psi[88,t,(s*6)] <- phiEl[4]*psiAbEl[s]
    } #s
    
    for (ew in 50:53){
    psi[85,t,ew] <- 0}
    psi[85,t,54] <- phiEl[1]*psi2El[9]*p9[3,t]*g9[t, 3]
    psi[85,t,55] <- phiEl[1]*psi2El[9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[85,t,56] <- phiEl[1]*psi2El[9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[85,t,57] <- phiEl[1]*psi2El[9]*(1-p9[3,t])
    for (ew in c(58:61, 63:66, 68:71)){
    psi[85,t,ew] <- 0}
    psi[85,t,72] <- phiEl[1]*psi2El[10]*p10[3,t]*g10[t, 3]
    psi[85,t,73] <- phiEl[1]*psi2El[10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[85,t,74] <- phiEl[1]*psi2El[10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[85,t,75] <- phiEl[1]*psi2El[10]*(1-p10[3,t])
    for (er in c(76:79,81:85,87:88)){
    psi[85,t,er] <- 0  }
    psi[85,t,86] <- phiEl[1]*psi2El[11]
    
    for (ed in 50:57){
    psi[86,t,ed] <- 0}
    psi[86,t,58] <- phiEl[2]*psi3El[9]*p9[4,t]*g9[t, 4]
    psi[86,t,59] <- phiEl[2]*psi3El[9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[86,t,60] <- phiEl[2]*psi3El[9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[86,t,61] <- phiEl[2]*psi3El[9]*(1-p9[4,t])
    for (ed in c(63:66,68:75)){
    psi[86,t,ed] <- 0}
    psi[86,t,76] <- phiEl[2]*psi3El[10]*p10[4,t]*g10[t, 4]
    psi[86,t,77] <- phiEl[2]*psi3El[10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[86,t,78] <- phiEl[2]*psi3El[10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[86,t,79] <- phiEl[2]*psi3El[10]*(1-p10[4,t])
    for (ef in 81:86){
    psi[86,t,ef] <- 0 }
    psi[86,t,87] <- phiEl[2]*psi3El[11]
    psi[86,t,88] <- 0
    
    for (es in 50:61){
    psi[87,t,es] <- 0}
    psi[87,t,63] <- phiEl[3]*psi4El[9]*p9[6,t]*g9[t, 6]
    psi[87,t,64] <- phiEl[3]*psi4El[9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[87,t,65] <- phiEl[3]*psi4El[9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[87,t,66] <- phiEl[3]*psi4El[9]*(1-p9[6,t])
    for (es in 68:79){
    psi[87,t,es] <- 0}
    psi[87,t,81] <- phiEl[3]*psi4El[10]*p10[6,t]*g10[t, 6]
    psi[87,t,82] <- phiEl[3]*psi4El[10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[87,t,83] <- phiEl[3]*psi4El[10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[87,t,84] <- phiEl[3]*psi4El[10]*(1-p10[6,t])
    psi[87,t,85] <- 0
    psi[87,t,86] <- 0
    psi[87,t,87] <- 0
    psi[87,t,88] <- phiEl[3]*psi4El[11]
    
    for (es in 50:61){
    psi[88,t,es] <- 0}
    psi[88,t,63] <- phiEl[4]*psiAbEl[9]*p9[6,t]*g9[t, 6]
    psi[88,t,64] <- phiEl[4]*psiAbEl[9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[88,t,65] <- phiEl[4]*psiAbEl[9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[88,t,66] <- phiEl[4]*psiAbEl[9]*(1-p9[6,t])
    for (es in 68:79){
    psi[88,t,es] <- 0}
    psi[88,t,81] <- phiEl[4]*psiAbEl[10]*p10[6,t]*g10[t, 6]
    psi[88,t,82] <- phiEl[4]*psiAbEl[10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[88,t,83] <- phiEl[4]*psiAbEl[10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[88,t,84] <- phiEl[4]*psiAbEl[10]*(1-p10[6,t])
    psi[88,t,85] <- 0
    psi[88,t,86] <- 0
    psi[88,t,87] <- 0
    psi[88,t,88] <- phiEl[4]*psiAbEl[11]
    
    
    # Define probabilities of O(t)
    for (j in 1:8){
    for (cc in c(2:4,6)){
    po[((j-1)*6)+cc,t] <- p[j,cc,t] }
    po[1+(6*(j-1)),t] <- 0 #State that can only occur on first capture
    po[5+(6*(j-1)),t] <- 0 #State that can only occur on first capture
    }
    for (jj in c(50:52, 54:56, 58:60, 63:65, 68:70, 72:74, 76:78, 81:83)){
    po[jj,t] <- 1
    }
    po[49,t] <- 0 #State that can only occur on first capture
    po[62,t] <- 0 #State that can only occur on first capture
    po[53,t] <- 0 #Unobservable state
    po[57,t] <- 0 #Unobservable state
    po[61,t] <- 0 #Unobservable state
    po[66,t] <- 0 #Unobservable state
    po[67,t] <- 0 #State that can only occur on first capture
    po[80,t] <- 0 #State that can only occur on first capture
    po[71,t] <- 0 #Unobservable state
    po[75,t] <- 0 #Unobservable state
    po[79,t] <- 0 #Unobservable state
    po[84,t] <- 0 #Unobservable state
    for (j in 85:88){
    po[j,t] <- pEl[j-84]
    }
    # Calculate probability of non-encounter (dq) and reshape the array for the encounter probabilities      
    for (s in 1:ns){
    dp[s,t,s] <- po[s,t]
    dq[s,t,s] <- 1-po[s,t]
    } # s
    
    for (s in 1:(ns-1)){
    for (m in (s+1):ns){
    dp[s,t,m] <- 0
    dq[s,t,m] <- 0
    } # s
    } # m
    for (s in 2:ns){
    for (m in 1:(s-1)){
    dp[s,t,m] <- 0
    dq[s,t,m] <- 0
    } # s
    } # m
    } # t
    
    # Define the multinomial likelihood
    for (t in 1:((n.occasions-1)*ns)){
    marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
    }
    
    # Define the cell probabilities of the multistate m-array   
    # Define matrix U: product of probabilities of state-transition and non-encounter (this is just done because there is no product function for matrix multiplication in JAGS)
    for (t in 1:(n.occasions-2)){
    U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
    for (j in (t+1):(n.occasions-1)){
    U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,t,] %*% dq[,t,]
    }
    }
    U[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ones
    
    # Diagonal
    for (t in 1:(n.occasions-2)){
    pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*% dp[,t,]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% psi[,j,] %*% dp[,j,]
    }
    }
    pr[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- psi[,n.occasions-1,] %*% dp[,n.occasions-1,]
    
    # Below main diagonal
    for (t in 2:(n.occasions-1)){
    for (j in 1:(t-1)){
    pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
    } #j
    } #t
    
    # Last column: probability of non-recapture
    for (t in 1:((n.occasions-1)*ns)){
    pr[t,(n.occasions*ns-(ns-1))] <- 1-sum(pr[t,1:((n.occasions-1)*ns)])
    } #t
    }    ",fill = TRUE)
sink()

# Bundle data
n.unobs <- 8 #Number of unobservable states
n.notobs <- 16 #Number of observable states that dont appear anywhere in capture histories
ns <- length(unique(as.numeric(mat.enc2))) - 1 + n.unobs + n.notobs # calculate the number of states

jags.data <- list(marr = ms.arr, n.occasions = ncol(mat.enc), rel = rowSums(ms.arr), ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), ind=ind, p.all.years=p.all.years)

inits <- function(){list(a.del=rep(0.5, length.out=3),  phiEl=rep(0.5, length.out=4), 
                         mean.phi = 0.1, sigma.phi = runif(1, 0, 5),
                         mean.p = 0.1, sigma.p = runif(1, 0, 5), A.delp=rnorm(3),
                         t.g9=rnorm(4), a.g.3 = 0.5,a.g.4 = 0.5,a.g.6 = 0.5, t.g10 = 0.5,
                         pEl=c(0.1,0.1,0.1,0.1),
                         theta2 = rnorm(1), theta3 = rnorm(1), theta4 = rnorm(1), thetaAa = rnorm(1), thetaAb = rnorm(1), repuls=rnorm(10), attract=rnorm(11),
                         mean.psiEl = 0.1)}  

# Parameters monitored
parameters <- c("phi.real", "sigma2.phi.real", "phiEl", "p.real", "sigma2.p.real", "p9.real","g9.real", "p10.real", "g10.real","pEl", "psi1", "psi2","psi3", "psiAa", "psiAb", "psi2El", "psiAbEl", "mean.psiEl")

# MCMC settings, very few iterations so any error shows quickly
ni <- 50 #000 
nt <- 1
nb <- 10 #000 
nc <- 3

# Call JAGS from R, using jagsUI syntax
starttime <- Sys.time()
m7 <- jags(jags.data, inits, parameters, "model7.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=FALSE)
endtime <- Sys.time()
options(max.print=999999) #Default is 99999
print(m7, digits = 3) #
print(m7$summary[1:400,], digits = 3)





########################
########################
#Model 8
########################
#Including a distance function for psi

#Same structure for phi and p as in model 7
#Same data as in model 7

# Analysis of the model
# Specify model in BUGS language
sink("model8.jags") 
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi[fro, stage, t]: survival probability at age/stage, having been at site fro, in period t
    # psi[fro, till]: movement probability from site fro to site till, only 10 of these as stage at next time-step is known and thereby constraints options to the 10 sites
    # p[...]: encounter probability. 
    # g: proportion of birds encountered (observed, captured or observed & captured) that where captured
    # -------------------------------------------------
    # States (S):
    # 1 alive at a1 s1
    # 2 alive at a2 s1
    # 3 alive at a3 s1
    # 4 alive at a4 s1
    # 5 alive at Ad_recent s1
    # 6 alive at Ad_previous s1
    # 7 alive at a1 s2
    # ...
    # 48 alive at Ad_previous s8
    # 49 alive & captured at a1 s9
    # 50 alive & captured at a2 s9
    # 51 alive & observed at a2 s9
    # 52 alive & observed & captured at a2 s9
    # 53 alive but not seen at a2
    # 54 alive & captured at a3 s10
    # ...
    # 85 alive Elsewhere at a2
    # 86 alive Elsewhere at a3
    # 87 alive Elsewhere at a4
    # 88 alive Elsewhere at Ad_previous
    # 89 dead
    #
    # Observations (O):
    # 1 seen at a1 s1
    # 2 seen at a2 s1
    # 3 seen at a3 s1
    # 4 seen at a4 s1
    # 5 seen at Ad s1
    # 6 seen at a1 s2
    # ...
    # 89/0 not seen
    # Unobservable states: 53,57,61,66, 71,75,79,84
    
    # -------------------------------------------------
    
    # Priors and constraints
    # Survival: random effects of time and island, additive effect of age
    for (fro in 1:10){
    for (stage in 1:6){ 
    for (t in 1:(n.occasions-1)){
    logit(phi[fro, stage, t]) <- mu.phi + eps.phi[fro, t] + del[stage] 
    }}}
    for (fro in 1:10){
    for (t in 1:(n.occasions-1)){
    eps.phi[fro, t] ~ dnorm(0, tau.phi)}}
    mean.phi ~ dunif(0, 1)                    # Prior for mean survival
    mu.phi <- log(mean.phi / (1-mean.phi))    # Logit transformation
    sigma.phi ~ dunif(0, 5)                   # Prior for standard deviation
    tau.phi <- pow(sigma.phi, -2)
    sigma2.phi <- pow(sigma.phi, 2)
    sigma2.phi.real <- sigma2.phi * pow(mean.phi, 2) * pow((1-mean.phi), 2) # Temporal variance on real scale
    
    del[1] <- 0                             # Priors for age/stage-effect
    for (stag in 2:4){ 
    del[stag] <- a.del[1]}
    a.del[1] ~ dnorm(0, 0.01)
    for (stag in 5:6){ 
    del[stag] <- a.del[stag-3]
    a.del[stag-3] ~ dnorm(0, 0.01)}
    for (fron in 1:10){  #Back-transform to the probability scale
    for (stage in 1:6){ 
    for (t in 1:(n.occasions-1)){
    phi.real[fron, stage, t] <- 1/(1+exp(-(mu.phi + eps.phi[fron,t] + del[stage])))    }}}
    for (st in 1:4){
    phiEl[st] ~ dunif(0, 1)}
    
    # Encounters at site 1-8: random effect of time and island, additive effect of age
    #Sites 1,3,8: 
    for (isl in p.all.years){ 
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[isl, stage, t]) <- mu.p + eps.p[isl, t] + delp[stage]     }}}
    for (isl in p.all.years){
    for (t in 1:(n.occasions-1)){
    eps.p[isl, t] ~ dnorm(0, tau.p)}}
    for (isln in p.all.years){  #Back-transform to the probability scale
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    p.real[isln, stage, t] <- 1/(1+exp(-(mu.p + eps.p[isln,t] + delp[stage])))    }}}
    #Site 2 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[2, stage, t]) <- mu.p + eps.p[2, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[2, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[2, stage, t] <- 1/(1+exp(-(mu.p + eps.p[2,t] + delp[stage])))    }}
    #Site 4 (should have p estimated for 2009-2011,2013-2017):
    for (stage in c(2:4,6)){ 
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    logit(p[4, stage, t]) <- mu.p + eps.p[4, t] + delp[stage]   }}
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    eps.p[4, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    p.real[4, stage, t] <- 1/(1+exp(-(mu.p + eps.p[4,t] + delp[stage])))    }}
    for (stage in c(2:4,6)){ 
    for (t in (n.occasions-6)){
    p[4, stage, t] <- 0   }}
    #Site 5 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[5, stage, t]) <- mu.p + eps.p[5, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[5, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[5, stage, t] <- 1/(1+exp(-(mu.p + eps.p[5,t] + delp[stage])))    }}
    #Site 6 (should have p estimated for 2011-2017):
    for (stage in c(2:4,6)){ 
    for (t in (n.occasions-7):(n.occasions-1)){
    logit(p[6, stage, t]) <- mu.p + eps.p[6, t] + delp[stage]   }}
    for (t in (n.occasions-7):(n.occasions-1)){
    eps.p[6, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in (n.occasions-7):(n.occasions-1)){
    p.real[6, stage, t] <- 1/(1+exp(-(mu.p + eps.p[6,t] + delp[stage])))    }}
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-8)){
    p[6, stage, t] <- 0   }}
    #Site 7 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[7, stage, t]) <- mu.p + eps.p[7, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[7, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[7, stage, t] <- 1/(1+exp(-(mu.p + eps.p[7,t] + delp[stage])))    }}
    
    mean.p ~ dunif(0, 1)                    # Prior for mean re-encounter
    mu.p <- log(mean.p / (1-mean.p))        # Logit transformation
    sigma.p ~ dunif(0, 5)                   # Prior for standard deviation
    tau.p <- pow(sigma.p, -2)
    sigma2.p <- pow(sigma.p, 2)
    sigma2.p.real <- sigma2.p * pow(mean.p, 2) * pow((1-mean.p), 2) # Temporal variance on real scale
    #delp[1] <- 0 #Re-encounters do not happen for age 1 [Shouldn't be necessary to define when using stage c(2:4,6) everywhere for p]
    #delp[5] <- 0 #Re-encounters do not happen for full-grown, just marked [Shouldn't be necessary to define when using stage c(2:4,6) everywhere for p]
    delp[2] <- 0 #Baseline age-class, which mu.p applies to
    for (stag in c(3,4)){ # Priors for age/stage-effect
    delp[stag] <- A.delp[stag-2]
    A.delp[stag-2] ~ dnorm(0, 0.01)}
    delp[6] <- A.delp[3]
    A.delp[3] ~ dnorm(0, 0.01)
    
    # Encounters at site 9-10:
    #Site 9 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p9[stage, t]) <- mu.p + eps.p9[t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p9[t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p9.real[stage, t] <- 1/(1+exp(-(mu.p + eps.p9[t] + delp[stage])))    }}
    #Site 10 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p10[stage, t]) <- mu.p + eps.p10[t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p10[t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p10.real[stage, t] <- 1/(1+exp(-(mu.p + eps.p10[t] + delp[stage])))    }}
    
    for (fd in 1:(n.occasions-5)){
    for (st in c(2:4,6)){
    g9[fd, st] <- 1  }}
    for (fd in (n.occasions-4):(n.occasions-1)){
    logit(g9[fd, 2]) <- t.g9[fd-n.occasions+5]
    logit(g9[fd, 3]) <- t.g9[fd-n.occasions+5] + a.g.3
    logit(g9[fd, 4]) <- t.g9[fd-n.occasions+5] + a.g.4
    logit(g9[fd, 6]) <- t.g9[fd-n.occasions+5] + a.g.6
    t.g9[fd-n.occasions+5] ~ dnorm(0, 0.01)
    g9.real[fd,2] <- 1/(1+exp(-(t.g9[fd-n.occasions+5]))) 
    g9.real[fd,3] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.3)))
    g9.real[fd,4] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.4)))
    g9.real[fd,6] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.6)))}
    for (fd in 1:(n.occasions-4)){ 
    for (st in c(2:4,6)){
    g10[fd, st] <- 1  }}
    logit(g10[n.occasions-3, 2]) <- t.g10
    logit(g10[n.occasions-3, 3]) <- t.g10 + a.g.3
    logit(g10[n.occasions-3, 4]) <- t.g10 + a.g.4
    logit(g10[n.occasions-3, 6]) <- t.g10 + a.g.6
    t.g10 ~ dnorm(0, 0.01)
    g10.real[n.occasions-3,2] <- 1/(1+exp(-(t.g10))) 
    g10.real[n.occasions-3,3] <- 1/(1+exp(-(t.g10 + a.g.3)))
    g10.real[n.occasions-3,4] <- 1/(1+exp(-(t.g10 + a.g.4)))
    g10.real[n.occasions-3,6] <- 1/(1+exp(-(t.g10 + a.g.6)))
    for (fd in (n.occasions-2):(n.occasions-1)){ 
    for (st in c(2:4,6)){
    g10[fd, st] <- 1  }}
    a.g.3 ~ dnorm(0, 0.01)
    a.g.4 ~ dnorm(0, 0.01)
    a.g.6 ~ dnorm(0, 0.01)
    
    # Encounters at site 11:
    for (st in 1:4){ 
    pEl[st] ~ dunif(0, 1) }
    
    # Transitions: multinomial logit
    # Normal priors on logit of all but one transition probs
    for (fro in 1:10){
    for (till in 1:(fro-1)){
    lpsi1[fro, till] <- repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi2[fro, till] <- theta2 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi3[fro, till] <- theta3 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi4[fro, till] <- theta4 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAa[fro, till] <- thetaAa + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAb[fro, till] <- thetaAb + repuls[fro] + attract[till] + B*dista[fro,till]}
    for (till in (fro+1):10){
    lpsi1[fro, till] <- repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi2[fro, till] <- theta2 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi3[fro, till] <- theta3 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi4[fro, till] <- theta4 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAa[fro, till] <- thetaAa + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAb[fro, till] <- thetaAb + repuls[fro] + attract[till] + B*dista[fro,till]} 
    lpsi1[fro, 11] <- repuls[fro] + attract[11]
    lpsi2[fro, 11] <- theta2 + repuls[fro] + attract[11]
    lpsi3[fro, 11] <- theta3 + repuls[fro] + attract[11]
    lpsi4[fro, 11] <- theta4 + repuls[fro] + attract[11]
    lpsiAa[fro, 11] <- thetaAa + repuls[fro] + attract[11]
    lpsiAb[fro, 11] <- thetaAb + repuls[fro] + attract[11]}
    
    for (fro in 1:10){
    for (till in 1:(fro-1)){
    psi1[fro, till] <- exp(lpsi1[fro, till]) / (1 + exp(lpsi1[fro,ind[fro,1]]) + exp(lpsi1[fro,ind[fro,2]]) + exp(lpsi1[fro,ind[fro,3]]) + exp(lpsi1[fro,ind[fro,4]]) + exp(lpsi1[fro,ind[fro,5]]) + exp(lpsi1[fro,ind[fro,6]]) + exp(lpsi1[fro,ind[fro,7]]) + exp(lpsi1[fro,ind[fro,8]]) + exp(lpsi1[fro,ind[fro,9]]) + exp(lpsi1[fro,ind[fro,10]]))
    psi2[fro, till] <- exp(lpsi2[fro, till]) / (1 + exp(lpsi2[fro,ind[fro,1]]) + exp(lpsi2[fro,ind[fro,2]]) + exp(lpsi2[fro,ind[fro,3]]) + exp(lpsi2[fro,ind[fro,4]]) + exp(lpsi2[fro,ind[fro,5]]) + exp(lpsi2[fro,ind[fro,6]]) + exp(lpsi2[fro,ind[fro,7]]) + exp(lpsi2[fro,ind[fro,8]]) + exp(lpsi2[fro,ind[fro,9]]) + exp(lpsi2[fro,ind[fro,10]]))
    psi3[fro, till] <- exp(lpsi3[fro, till]) / (1 + exp(lpsi3[fro,ind[fro,1]]) + exp(lpsi3[fro,ind[fro,2]]) + exp(lpsi3[fro,ind[fro,3]]) + exp(lpsi3[fro,ind[fro,4]]) + exp(lpsi3[fro,ind[fro,5]]) + exp(lpsi3[fro,ind[fro,6]]) + exp(lpsi3[fro,ind[fro,7]]) + exp(lpsi3[fro,ind[fro,8]]) + exp(lpsi3[fro,ind[fro,9]]) + exp(lpsi3[fro,ind[fro,10]]))
    psi4[fro, till] <- exp(lpsi4[fro, till]) / (1 + exp(lpsi4[fro,ind[fro,1]]) + exp(lpsi4[fro,ind[fro,2]]) + exp(lpsi4[fro,ind[fro,3]]) + exp(lpsi4[fro,ind[fro,4]]) + exp(lpsi4[fro,ind[fro,5]]) + exp(lpsi4[fro,ind[fro,6]]) + exp(lpsi4[fro,ind[fro,7]]) + exp(lpsi4[fro,ind[fro,8]]) + exp(lpsi4[fro,ind[fro,9]]) + exp(lpsi4[fro,ind[fro,10]]))
    psiAa[fro, till] <- exp(lpsiAa[fro, till]) / (1 + exp(lpsiAa[fro,ind[fro,1]]) + exp(lpsiAa[fro,ind[fro,2]]) + exp(lpsiAa[fro,ind[fro,3]]) + exp(lpsiAa[fro,ind[fro,4]]) + exp(lpsiAa[fro,ind[fro,5]]) + exp(lpsiAa[fro,ind[fro,6]]) + exp(lpsiAa[fro,ind[fro,7]]) + exp(lpsiAa[fro,ind[fro,8]]) + exp(lpsiAa[fro,ind[fro,9]]) + exp(lpsiAa[fro,ind[fro,10]]))
    psiAb[fro, till] <- exp(lpsiAb[fro, till]) / (1 + exp(lpsiAb[fro,ind[fro,1]]) + exp(lpsiAb[fro,ind[fro,2]]) + exp(lpsiAb[fro,ind[fro,3]]) + exp(lpsiAb[fro,ind[fro,4]]) + exp(lpsiAb[fro,ind[fro,5]]) + exp(lpsiAb[fro,ind[fro,6]]) + exp(lpsiAb[fro,ind[fro,7]]) + exp(lpsiAb[fro,ind[fro,8]]) + exp(lpsiAb[fro,ind[fro,9]]) + exp(lpsiAb[fro,ind[fro,10]]))  }
    for (till in (fro+1):11){
    psi1[fro, till] <- exp(lpsi1[fro, till]) / (1 + exp(lpsi1[fro,ind[fro,1]]) + exp(lpsi1[fro,ind[fro,2]]) + exp(lpsi1[fro,ind[fro,3]]) + exp(lpsi1[fro,ind[fro,4]]) + exp(lpsi1[fro,ind[fro,5]]) + exp(lpsi1[fro,ind[fro,6]]) + exp(lpsi1[fro,ind[fro,7]]) + exp(lpsi1[fro,ind[fro,8]]) + exp(lpsi1[fro,ind[fro,9]]) + exp(lpsi1[fro,ind[fro,10]]))
    psi2[fro, till] <- exp(lpsi2[fro, till]) / (1 + exp(lpsi2[fro,ind[fro,1]]) + exp(lpsi2[fro,ind[fro,2]]) + exp(lpsi2[fro,ind[fro,3]]) + exp(lpsi2[fro,ind[fro,4]]) + exp(lpsi2[fro,ind[fro,5]]) + exp(lpsi2[fro,ind[fro,6]]) + exp(lpsi2[fro,ind[fro,7]]) + exp(lpsi2[fro,ind[fro,8]]) + exp(lpsi2[fro,ind[fro,9]]) + exp(lpsi2[fro,ind[fro,10]]))
    psi3[fro, till] <- exp(lpsi3[fro, till]) / (1 + exp(lpsi3[fro,ind[fro,1]]) + exp(lpsi3[fro,ind[fro,2]]) + exp(lpsi3[fro,ind[fro,3]]) + exp(lpsi3[fro,ind[fro,4]]) + exp(lpsi3[fro,ind[fro,5]]) + exp(lpsi3[fro,ind[fro,6]]) + exp(lpsi3[fro,ind[fro,7]]) + exp(lpsi3[fro,ind[fro,8]]) + exp(lpsi3[fro,ind[fro,9]]) + exp(lpsi3[fro,ind[fro,10]]))
    psi4[fro, till] <- exp(lpsi4[fro, till]) / (1 + exp(lpsi4[fro,ind[fro,1]]) + exp(lpsi4[fro,ind[fro,2]]) + exp(lpsi4[fro,ind[fro,3]]) + exp(lpsi4[fro,ind[fro,4]]) + exp(lpsi4[fro,ind[fro,5]]) + exp(lpsi4[fro,ind[fro,6]]) + exp(lpsi4[fro,ind[fro,7]]) + exp(lpsi4[fro,ind[fro,8]]) + exp(lpsi4[fro,ind[fro,9]]) + exp(lpsi4[fro,ind[fro,10]]))
    psiAa[fro, till] <- exp(lpsiAa[fro, till]) / (1 + exp(lpsiAa[fro,ind[fro,1]]) + exp(lpsiAa[fro,ind[fro,2]]) + exp(lpsiAa[fro,ind[fro,3]]) + exp(lpsiAa[fro,ind[fro,4]]) + exp(lpsiAa[fro,ind[fro,5]]) + exp(lpsiAa[fro,ind[fro,6]]) + exp(lpsiAa[fro,ind[fro,7]]) + exp(lpsiAa[fro,ind[fro,8]]) + exp(lpsiAa[fro,ind[fro,9]]) + exp(lpsiAa[fro,ind[fro,10]]))
    psiAb[fro, till] <- exp(lpsiAb[fro, till]) / (1 + exp(lpsiAb[fro,ind[fro,1]]) + exp(lpsiAb[fro,ind[fro,2]]) + exp(lpsiAb[fro,ind[fro,3]]) + exp(lpsiAb[fro,ind[fro,4]]) + exp(lpsiAb[fro,ind[fro,5]]) + exp(lpsiAb[fro,ind[fro,6]]) + exp(lpsiAb[fro,ind[fro,7]]) + exp(lpsiAb[fro,ind[fro,8]]) + exp(lpsiAb[fro,ind[fro,9]]) + exp(lpsiAb[fro,ind[fro,10]]))  }
    #The case below is when fro=till, but rewritten to match JAGS syntax
    psi1[fro,fro] <- 1- sum(psi1[fro,ind[fro,]]) 
    psi2[fro,fro] <- 1- sum(psi2[fro,ind[fro,]]) 
    psi3[fro,fro] <- 1- sum(psi3[fro,ind[fro,]]) 
    psi4[fro,fro] <- 1- sum(psi4[fro,ind[fro,]]) 
    psiAa[fro,fro] <- 1- sum(psiAa[fro,ind[fro,]])   
    psiAb[fro,fro] <- 1- sum(psiAb[fro,ind[fro,]])   
    } 
    
    #Priors for transitions
    theta2 ~ dnorm(0, 0.01)               # Prior for age-effect in psi
    theta3 ~ dnorm(0, 0.01)  
    theta4 ~ dnorm(0, 0.01)  
    thetaAa ~ dnorm(0, 0.01)  
    thetaAb ~ dnorm(0, 0.01) 
    #mean.psi ~ dunif(0, 1)                # Prior for mean movement prob - not used as no transition makes sense to set as baseline
    #mu <- log(mean.psi / (1-mean.psi))    # Logit transformation
    for (w in 1:10){
    repuls[w] ~ dnorm(0, 0.01) 
    attract[w] ~ dnorm(0, 0.01) 
    }
    B ~ dnorm(0, 0.001)I(-10,10)
    attract[11] ~ dnorm(0, 0.01)
    
    #Transitions from Elsewhere
    for (till in 1:10){
    lpsi2El[till] <- muEl 
    lpsi3El[till] <- muEl + theta3
    lpsi4El[till] <- muEl + theta4
    lpsiAbEl[till] <- muEl + thetaAb
    }
    for (till in 1:10){
    psi2El[till] <- exp(lpsi2El[till]) / (1 + exp(lpsi2El[1]) + exp(lpsi2El[2]) + exp(lpsi2El[3]) + exp(lpsi2El[4]) + exp(lpsi2El[5]) + exp(lpsi2El[6]) + exp(lpsi2El[7]) + exp(lpsi2El[8]) + exp(lpsi2El[9]) + exp(lpsi2El[10]))
    psi3El[till] <- exp(lpsi3El[till]) / (1 + exp(lpsi3El[1]) + exp(lpsi3El[2]) + exp(lpsi3El[3]) + exp(lpsi3El[4]) + exp(lpsi3El[5]) + exp(lpsi3El[6]) + exp(lpsi3El[7]) + exp(lpsi3El[8]) + exp(lpsi3El[9]) + exp(lpsi3El[10]))
    psi4El[till] <- exp(lpsi4El[till]) / (1 + exp(lpsi4El[1]) + exp(lpsi4El[2]) + exp(lpsi4El[3]) + exp(lpsi4El[4]) + exp(lpsi4El[5]) + exp(lpsi4El[6]) + exp(lpsi4El[7]) + exp(lpsi4El[8]) + exp(lpsi4El[9]) + exp(lpsi4El[10]))
    psiAbEl[till] <- exp(lpsiAbEl[till]) / (1 + exp(lpsiAbEl[1]) + exp(lpsiAbEl[2]) + exp(lpsiAbEl[3]) + exp(lpsiAbEl[4]) + exp(lpsiAbEl[5]) + exp(lpsiAbEl[6]) + exp(lpsiAbEl[7]) + exp(lpsiAbEl[8]) + exp(lpsiAbEl[9]) + exp(lpsiAbEl[10]))
    }
    # Calculate the last transition probability
    psi2El[11] <- 1- sum(psi2El[1:10]) 
    psi3El[11] <- 1- sum(psi3El[1:10]) 
    psi4El[11] <- 1- sum(psi4El[1:10]) 
    psiAbEl[11] <- 1- sum(psiAbEl[1:10])   
    #Prior for transitions from Elsewhere    
    mean.psiEl ~ dunif(0, 1)                   # Prior for mean movement prob from Elsewhere
    muEl <- log(mean.psiEl / (1-mean.psiEl))     # Logit transformation
    
    
    # Define state-transition and observation matrices 	
    # Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
    for (d in 1:88){
    for (m in 1:8){
    psi[d,t,1+(6*(m-1))] <- 0
    psi[d,t,5+(6*(m-1))] <- 0
    } #m
    psi[d,t,49] <- 0
    psi[d,t,62] <- 0
    psi[d,t,67] <- 0
    psi[d,t,80] <- 0
    } #d
    
    for (sf in 1:8){  
    for (s in 1:8){
    psi[(6*(sf-1))+1,t,(s*6-4)] <- phi[sf,1,t]*psi1[sf,s]
    psi[(6*(sf-1))+1,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+1,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+1,t,(s*6)] <- 0
    psi[(6*(sf-1))+2,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+2,t,(s*6-3)] <- phi[sf,2,t]*psi2[sf,s]
    psi[(6*(sf-1))+2,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+2,t,(s*6)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-2)] <- phi[sf,3,t]*psi3[sf,s]
    psi[(6*(sf-1))+3,t,(s*6)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+4,t,(s*6)] <- phi[sf,4,t]*psi4[sf,s]
    psi[(6*(sf-1))+5,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+5,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+5,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+5,t,(s*6)] <- phi[sf,5,t]*psiAa[sf,s]
    psi[(6*(sf-1))+6,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+6,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+6,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+6,t,(s*6)] <- phi[sf,6,t]*psiAb[sf,s]
    } #s
    
    psi[(6*(sf-1))+1,t,50] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*g9[t, 2]
    psi[(6*(sf-1))+1,t,51] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[(6*(sf-1))+1,t,52] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[(6*(sf-1))+1,t,53] <- phi[sf,1,t]*psi1[sf,9]*(1-p9[2,t])
    for (ex in c(54:61,63:66)){
    psi[(6*(sf-1))+1,t,ex] <- 0  }
    psi[(6*(sf-1))+1,t,68] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*g10[t, 2]
    psi[(6*(sf-1))+1,t,69] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[(6*(sf-1))+1,t,70] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[(6*(sf-1))+1,t,71] <- phi[sf,1,t]*psi1[sf,10]*(1-p10[2,t])
    for (ex in c(72:79,81:84)){
    psi[(6*(sf-1))+1,t,ex] <- 0  }
    psi[(6*(sf-1))+1,t,85] <- phi[sf,1,t]*psi1[sf,11]
    psi[(6*(sf-1))+1,t,86] <- 0
    psi[(6*(sf-1))+1,t,87] <- 0
    psi[(6*(sf-1))+1,t,88] <- 0
    
    for (ew in 50:53){
    psi[(6*(sf-1))+2,t,ew] <- 0}
    psi[(6*(sf-1))+2,t,54] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*g9[t, 3]
    psi[(6*(sf-1))+2,t,55] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[(6*(sf-1))+2,t,56] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[(6*(sf-1))+2,t,57] <- phi[sf,2,t]*psi2[sf,9]*(1-p9[3,t])
    for (er in c(58:61,63:66)){
    psi[(6*(sf-1))+2,t,er] <- 0  }
    for (ew in 68:71){
    psi[(6*(sf-1))+2,t,ew] <- 0}
    psi[(6*(sf-1))+2,t,72] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*g10[t, 3]
    psi[(6*(sf-1))+2,t,73] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[(6*(sf-1))+2,t,74] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[(6*(sf-1))+2,t,75] <- phi[sf,2,t]*psi2[sf,10]*(1-p10[3,t])
    for (er in c(76:79,81:84)){
    psi[(6*(sf-1))+2,t,er] <- 0  }
    psi[(6*(sf-1))+2,t,85] <- 0
    psi[(6*(sf-1))+2,t,86] <- phi[sf,2,t]*psi2[sf,11]
    psi[(6*(sf-1))+2,t,87] <- 0
    psi[(6*(sf-1))+2,t,88] <- 0
    
    for (ed in 50:57){
    psi[(6*(sf-1))+3,t,ed] <- 0}
    psi[(6*(sf-1))+3,t,58] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*g9[t, 4]
    psi[(6*(sf-1))+3,t,59] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[(6*(sf-1))+3,t,60] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[(6*(sf-1))+3,t,61] <- phi[sf,3,t]*psi3[sf,9]*(1-p9[4,t])
    for (ef in 63:66){
    psi[(6*(sf-1))+3,t,ef] <- 0 }
    for (ed in 68:75){
    psi[(6*(sf-1))+3,t,ed] <- 0}
    psi[(6*(sf-1))+3,t,76] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*g10[t, 4]
    psi[(6*(sf-1))+3,t,77] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[(6*(sf-1))+3,t,78] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[(6*(sf-1))+3,t,79] <- phi[sf,3,t]*psi3[sf,10]*(1-p10[4,t])
    for (ef in 81:84){
    psi[(6*(sf-1))+3,t,ef] <- 0 }
    psi[(6*(sf-1))+3,t,85] <- 0
    psi[(6*(sf-1))+3,t,86] <- 0
    psi[(6*(sf-1))+3,t,87] <- phi[sf,3,t]*psi3[sf,11]
    psi[(6*(sf-1))+3,t,88] <- 0
    
    for (es in 50:61){
    psi[(6*(sf-1))+4,t,es] <- 0}
    psi[(6*(sf-1))+4,t,63] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+4,t,64] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+4,t,65] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+4,t,66] <- phi[sf,4,t]*psi4[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+4,t,es] <- 0}
    psi[(6*(sf-1))+4,t,81] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+4,t,82] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+4,t,83] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+4,t,84] <- phi[sf,4,t]*psi4[sf,10]*(1-p10[6,t])
    psi[(6*(sf-1))+4,t,85] <- 0
    psi[(6*(sf-1))+4,t,86] <- 0
    psi[(6*(sf-1))+4,t,87] <- 0
    psi[(6*(sf-1))+4,t,88] <- phi[sf,4,t]*psi4[sf,11]
    
    for (es in 50:61){
    psi[(6*(sf-1))+5,t,es] <- 0}
    psi[(6*(sf-1))+5,t,63] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+5,t,64] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+5,t,65] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+5,t,66] <- phi[sf,5,t]*psiAa[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+5,t,es] <- 0}
    psi[(6*(sf-1))+5,t,81] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+5,t,82] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+5,t,83] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+5,t,84] <- phi[sf,5,t]*psiAa[sf,10]*(1-p10[6,t])
    psi[(6*(sf-1))+5,t,85] <- 0
    psi[(6*(sf-1))+5,t,86] <- 0
    psi[(6*(sf-1))+5,t,87] <- 0
    psi[(6*(sf-1))+5,t,88] <- phi[sf,5,t]*psiAa[sf,11]
    
    for (es in 50:61){
    psi[(6*(sf-1))+6,t,es] <- 0}
    psi[(6*(sf-1))+6,t,63] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+6,t,64] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+6,t,65] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+6,t,66] <- phi[sf,6,t]*psiAb[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+6,t,es] <- 0}
    psi[(6*(sf-1))+6,t,81] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+6,t,82] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+6,t,83] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+6,t,84] <- phi[sf,6,t]*psiAb[sf,10]*(1-p10[6,t])
    psi[(6*(sf-1))+6,t,85] <- 0
    psi[(6*(sf-1))+6,t,86] <- 0
    psi[(6*(sf-1))+6,t,87] <- 0
    psi[(6*(sf-1))+6,t,88] <- phi[sf,6,t]*psiAb[sf,11]
    } #sf
    
    
    for (s in 1:8){
    psi[49,t,(s*6-4)] <- phi[9,1,t]*psi1[9,s]
    psi[49,t,(s*6-3)] <- 0
    psi[49,t,(s*6-2)] <- 0
    psi[49,t,(s*6)] <- 0
    for (i in 1:4){
    psi[49+i,t,(s*6-4)] <- 0
    psi[49+i,t,(s*6-3)] <- phi[9,2,t]*psi2[9,s]
    psi[49+i,t,(s*6-2)] <- 0
    psi[49+i,t,(s*6)] <- 0
    psi[53+i,t,(s*6-4)] <- 0
    psi[53+i,t,(s*6-3)] <- 0
    psi[53+i,t,(s*6-2)] <- phi[9,3,t]*psi3[9,s]
    psi[53+i,t,(s*6)] <- 0
    psi[57+i,t,(s*6-4)] <- 0
    psi[57+i,t,(s*6-3)] <- 0
    psi[57+i,t,(s*6-2)] <- 0
    psi[57+i,t,(s*6)] <- phi[9,4,t]*psi4[9,s] 
    psi[62+i,t,(s*6-4)] <- 0
    psi[62+i,t,(s*6-3)] <- 0
    psi[62+i,t,(s*6-2)] <- 0
    psi[62+i,t,(s*6)] <- phi[9,6,t]*psiAb[9,s]
    } #i
    psi[62,t,(s*6-4)] <- 0
    psi[62,t,(s*6-3)] <- 0
    psi[62,t,(s*6-2)] <- 0
    psi[62,t,(s*6)] <- phi[9,5,t]*psiAa[9,s]
    } #s
    
    for (s in 1:8){
    psi[67,t,(s*6-4)] <- phi[10,1,t]*psi1[10,s]
    psi[67,t,(s*6-3)] <- 0
    psi[67,t,(s*6-2)] <- 0
    psi[67,t,(s*6)] <- 0
    for (i in 1:4){
    psi[67+i,t,(s*6-4)] <- 0
    psi[67+i,t,(s*6-3)] <- phi[10,2,t]*psi2[10,s]
    psi[67+i,t,(s*6-2)] <- 0
    psi[67+i,t,(s*6)] <- 0
    psi[71+i,t,(s*6-4)] <- 0
    psi[71+i,t,(s*6-3)] <- 0
    psi[71+i,t,(s*6-2)] <- phi[10,3,t]*psi3[10,s]
    psi[71+i,t,(s*6)] <- 0
    psi[75+i,t,(s*6-4)] <- 0
    psi[75+i,t,(s*6-3)] <- 0
    psi[75+i,t,(s*6-2)] <- 0
    psi[75+i,t,(s*6)] <- phi[10,4,t]*psi4[10,s] 
    psi[80+i,t,(s*6-4)] <- 0
    psi[80+i,t,(s*6-3)] <- 0
    psi[80+i,t,(s*6-2)] <- 0
    psi[80+i,t,(s*6)] <- phi[10,6,t]*psiAb[10,s]
    } #i
    psi[80,t,(s*6-4)] <- 0
    psi[80,t,(s*6-3)] <- 0
    psi[80,t,(s*6-2)] <- 0
    psi[80,t,(s*6)] <- phi[10,5,t]*psiAa[10,s]
    } #s
    
    psi[49,t,50] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*g9[t, 2]
    psi[49,t,51] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[49,t,52] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[49,t,53] <- phi[9,1,t]*psi1[9,9]*(1-p9[2,t])
    for (ex in c(54:61,63:66)){
    psi[49,t,ex] <- 0  }
    
    for (i in 1:4){
    for (ew in 50:53){
    psi[49+i,t,ew] <- 0  }
    psi[49+i,t,54] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*g9[t, 3]
    psi[49+i,t,55] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[49+i,t,56] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[49+i,t,57] <- phi[9,2,t]*psi2[9,9]*(1-p9[3,t])
    for (er in c(58:61,63:66)){
    psi[49+i,t,er] <- 0  }
    for (ed in 50:57){
    psi[53+i,t,ed] <- 0}
    psi[53+i,t,58] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*g9[t, 4]
    psi[53+i,t,59] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[53+i,t,60] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[53+i,t,61] <- phi[9,3,t]*psi3[9,9]*(1-p9[4,t])
    for (ef in 63:66){
    psi[53+i,t,ef] <- 0 }
    for (es in 50:61){
    psi[57+i,t,es] <- 0}
    psi[57+i,t,63] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*g9[t, 6]
    psi[57+i,t,64] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[57+i,t,65] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[57+i,t,66] <- phi[9,4,t]*psi4[9,9]*(1-p9[6,t])
    for (es in 50:61){
    psi[62+i,t,es] <- 0}
    psi[62+i,t,63] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*g9[t, 6]
    psi[62+i,t,64] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62+i,t,65] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62+i,t,66] <- phi[9,6,t]*psiAb[9,9]*(1-p9[6,t])
    } #i
    
    for (es in 50:61){
    psi[62,t,es] <- 0}
    psi[62,t,63] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*g9[t, 6]
    psi[62,t,64] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62,t,65] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62,t,66] <- phi[9,5,t]*psiAa[9,9]*(1-p9[6,t])
    
    
    psi[49,t,50+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*g10[t, 2]
    psi[49,t,51+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[49,t,52+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[49,t,53+18] <- phi[9,1,t]*psi1[9,10]*(1-p10[2,t])
    for (ex in c(54:61,63:66)){
    psi[49,t,ex+18] <- 0  }
    
    for (i in 1:4){
    for (ew in 50:53){
    psi[49+i,t,ew+18] <- 0  }
    psi[49+i,t,54+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*g10[t, 3]
    psi[49+i,t,55+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[49+i,t,56+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[49+i,t,57+18] <- phi[9,2,t]*psi2[9,10]*(1-p10[3,t])
    for (er in c(58:61,63:66)){
    psi[49+i,t,er+18] <- 0  }
    for (ed in 50:57){
    psi[53+i,t,ed+18] <- 0}
    psi[53+i,t,58+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*g10[t, 4]
    psi[53+i,t,59+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[53+i,t,60+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[53+i,t,61+18] <- phi[9,3,t]*psi3[9,10]*(1-p10[4,t])
    for (ef in 63:66){
    psi[53+i,t,ef+18] <- 0 }
    for (es in 50:61){
    psi[57+i,t,es+18] <- 0}
    psi[57+i,t,63+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*g10[t, 6]
    psi[57+i,t,64+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[57+i,t,65+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[57+i,t,66+18] <- phi[9,4,t]*psi4[9,10]*(1-p10[6,t])
    for (es in 50:61){
    psi[62+i,t,es+18] <- 0}
    psi[62+i,t,63+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*g10[t, 6]
    psi[62+i,t,64+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62+i,t,65+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62+i,t,66+18] <- phi[9,6,t]*psiAb[9,10]*(1-p10[6,t])
    } #i
    
    for (es in 50:61){
    psi[62,t,es+18] <- 0}
    psi[62,t,63+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*g10[t, 6]
    psi[62,t,64+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62,t,65+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62,t,66+18] <- phi[9,5,t]*psiAa[9,10]*(1-p10[6,t])
    
    
    psi[49,t,(89-4)] <- phi[9,1,t]*psi1[9,11]
    psi[49,t,(89-3)] <- 0
    psi[49,t,(89-2)] <- 0
    psi[49,t,(88)] <- 0
    for (i in 1:4){
    psi[49+i,t,(89-4)] <- 0
    psi[49+i,t,(89-3)] <- phi[9,2,t]*psi2[9,11]
    psi[49+i,t,(89-2)] <- 0
    psi[49+i,t,(88)] <- 0
    psi[53+i,t,(89-4)] <- 0
    psi[53+i,t,(89-3)] <- 0
    psi[53+i,t,(89-2)] <- phi[9,3,t]*psi3[9,11]
    psi[53+i,t,(88)] <- 0
    psi[57+i,t,(89-4)] <- 0
    psi[57+i,t,(89-3)] <- 0
    psi[57+i,t,(89-2)] <- 0
    psi[57+i,t,(88)] <- phi[9,4,t]*psi4[9,11] 
    psi[62+i,t,(89-4)] <- 0
    psi[62+i,t,(89-3)] <- 0
    psi[62+i,t,(89-2)] <- 0
    psi[62+i,t,(88)] <- phi[9,6,t]*psiAb[9,11]
    } #i
    psi[62,t,(89-4)] <- 0
    psi[62,t,(89-3)] <- 0
    psi[62,t,(89-2)] <- 0
    psi[62,t,(88)] <- phi[9,5,t]*psiAa[9,11]
    
    
    psi[67,t,50] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*g9[t, 2]
    psi[67,t,51] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[67,t,52] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[67,t,53] <- phi[10,1,t]*psi1[10,9]*(1-p9[2,t])
    for (ex in c(72:79,81:84)){
    psi[67,t,ex-18] <- 0  }
    
    for (i in 1:4){
    for (ew in 68:71){
    psi[67+i,t,ew-18] <- 0  }
    psi[67+i,t,72-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*g9[t, 3]
    psi[67+i,t,73-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[67+i,t,74-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[67+i,t,75-18] <- phi[10,2,t]*psi2[10,9]*(1-p9[3,t])
    for (er in c(76:79,81:84)){
    psi[67+i,t,er-18] <- 0  }
    for (ed in 68:75){
    psi[71+i,t,ed-18] <- 0}
    psi[71+i,t,76-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*g9[t, 4]
    psi[71+i,t,77-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[71+i,t,78-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[71+i,t,79-18] <- phi[10,3,t]*psi3[10,9]*(1-p9[4,t])
    for (ef in 81:84){
    psi[71+i,t,ef-18] <- 0 }
    for (es in 68:79){
    psi[75+i,t,es-18] <- 0}
    psi[75+i,t,81-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*g9[t, 6]
    psi[75+i,t,82-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[75+i,t,83-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[75+i,t,84-18] <- phi[10,4,t]*psi4[10,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[80+i,t,es-18] <- 0}
    psi[80+i,t,81-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*g9[t, 6]
    psi[80+i,t,82-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80+i,t,83-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80+i,t,84-18] <- phi[10,6,t]*psiAb[10,9]*(1-p9[6,t])
    } #i
    
    for (es in 68:79){
    psi[80,t,es-18] <- 0}
    psi[80,t,81-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*g9[t, 6]
    psi[80,t,82-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80,t,83-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80,t,84-18] <- phi[10,5,t]*psiAa[10,9]*(1-p9[6,t])
    
    
    psi[67,t,68] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*g10[t, 2]
    psi[67,t,69] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[67,t,70] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[67,t,71] <- phi[10,1,t]*psi1[10,10]*(1-p10[2,t])
    for (ex in c(72:79,81:84)){
    psi[67,t,ex] <- 0  }
    
    for (i in 1:4){
    for (ew in 68:71){
    psi[67+i,t,ew] <- 0  }
    psi[67+i,t,72] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*g10[t, 3]
    psi[67+i,t,73] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[67+i,t,74] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[67+i,t,75] <- phi[10,2,t]*psi2[10,10]*(1-p10[3,t])
    for (er in c(76:79,81:84)){
    psi[67+i,t,er] <- 0  }
    for (ed in 68:75){
    psi[71+i,t,ed] <- 0}
    psi[71+i,t,76] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*g10[t, 4]
    psi[71+i,t,77] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[71+i,t,78] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[71+i,t,79] <- phi[10,3,t]*psi3[10,10]*(1-p10[4,t])
    for (ef in 81:84){
    psi[71+i,t,ef] <- 0 }
    for (es in 68:79){
    psi[75+i,t,es] <- 0}
    psi[75+i,t,81] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*g10[t, 6]
    psi[75+i,t,82] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[75+i,t,83] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[75+i,t,84] <- phi[10,4,t]*psi4[10,10]*(1-p10[6,t])
    for (es in 68:79){
    psi[80+i,t,es] <- 0}
    psi[80+i,t,81] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*g10[t, 6]
    psi[80+i,t,82] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80+i,t,83] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80+i,t,84] <- phi[10,6,t]*psiAb[10,10]*(1-p10[6,t])
    } #i
    
    for (es in 68:79){
    psi[80,t,es] <- 0}
    psi[80,t,81] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*g10[t, 6]
    psi[80,t,82] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80,t,83] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80,t,84] <- phi[10,5,t]*psiAa[10,10]*(1-p10[6,t])
    
    
    psi[67,t,(89-4)] <- phi[10,1,t]*psi1[10,11]
    psi[67,t,(89-3)] <- 0
    psi[67,t,(89-2)] <- 0
    psi[67,t,(88)] <- 0
    for (i in 1:4){
    psi[67+i,t,(89-4)] <- 0
    psi[67+i,t,(89-3)] <- phi[10,2,t]*psi2[10,11]
    psi[67+i,t,(89-2)] <- 0
    psi[67+i,t,(88)] <- 0
    psi[71+i,t,(89-4)] <- 0
    psi[71+i,t,(89-3)] <- 0
    psi[71+i,t,(89-2)] <- phi[10,3,t]*psi3[10,11]
    psi[71+i,t,(88)] <- 0
    psi[75+i,t,(89-4)] <- 0
    psi[75+i,t,(89-3)] <- 0
    psi[75+i,t,(89-2)] <- 0
    psi[75+i,t,(88)] <- phi[10,4,t]*psi4[10,11] 
    psi[80+i,t,(89-4)] <- 0
    psi[80+i,t,(89-3)] <- 0
    psi[80+i,t,(89-2)] <- 0
    psi[80+i,t,(88)] <- phi[10,6,t]*psiAb[10,11]
    } #i
    psi[80,t,(89-4)] <- 0
    psi[80,t,(89-3)] <- 0
    psi[80,t,(89-2)] <- 0
    psi[80,t,(88)] <- phi[10,5,t]*psiAa[10,11]
    
    for (s in 1:8){
    psi[85,t,(s*6-4)] <- 0
    psi[85,t,(s*6-3)] <- phiEl[1]*psi2El[s]
    psi[85,t,(s*6-2)] <- 0
    psi[85,t,(s*6)] <- 0
    psi[86,t,(s*6-4)] <- 0
    psi[86,t,(s*6-3)] <- 0
    psi[86,t,(s*6-2)] <- phiEl[2]*psi3El[s]
    psi[86,t,(s*6)] <- 0
    psi[87,t,(s*6-4)] <- 0
    psi[87,t,(s*6-3)] <- 0
    psi[87,t,(s*6-2)] <- 0
    psi[87,t,(s*6)] <- phiEl[3]*psi4El[s]
    psi[88,t,(s*6-4)] <- 0
    psi[88,t,(s*6-3)] <- 0
    psi[88,t,(s*6-2)] <- 0
    psi[88,t,(s*6)] <- phiEl[4]*psiAbEl[s]
    } #s
    
    for (ew in 50:53){
    psi[85,t,ew] <- 0}
    psi[85,t,54] <- phiEl[1]*psi2El[9]*p9[3,t]*g9[t, 3]
    psi[85,t,55] <- phiEl[1]*psi2El[9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[85,t,56] <- phiEl[1]*psi2El[9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[85,t,57] <- phiEl[1]*psi2El[9]*(1-p9[3,t])
    for (ew in c(58:61, 63:66, 68:71)){
    psi[85,t,ew] <- 0}
    psi[85,t,72] <- phiEl[1]*psi2El[10]*p10[3,t]*g10[t, 3]
    psi[85,t,73] <- phiEl[1]*psi2El[10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[85,t,74] <- phiEl[1]*psi2El[10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[85,t,75] <- phiEl[1]*psi2El[10]*(1-p10[3,t])
    for (er in c(76:79,81:85,87:88)){
    psi[85,t,er] <- 0  }
    psi[85,t,86] <- phiEl[1]*psi2El[11]
    
    for (ed in 50:57){
    psi[86,t,ed] <- 0}
    psi[86,t,58] <- phiEl[2]*psi3El[9]*p9[4,t]*g9[t, 4]
    psi[86,t,59] <- phiEl[2]*psi3El[9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[86,t,60] <- phiEl[2]*psi3El[9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[86,t,61] <- phiEl[2]*psi3El[9]*(1-p9[4,t])
    for (ed in c(63:66,68:75)){
    psi[86,t,ed] <- 0}
    psi[86,t,76] <- phiEl[2]*psi3El[10]*p10[4,t]*g10[t, 4]
    psi[86,t,77] <- phiEl[2]*psi3El[10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[86,t,78] <- phiEl[2]*psi3El[10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[86,t,79] <- phiEl[2]*psi3El[10]*(1-p10[4,t])
    for (ef in 81:86){
    psi[86,t,ef] <- 0 }
    psi[86,t,87] <- phiEl[2]*psi3El[11]
    psi[86,t,88] <- 0
    
    for (es in 50:61){
    psi[87,t,es] <- 0}
    psi[87,t,63] <- phiEl[3]*psi4El[9]*p9[6,t]*g9[t, 6]
    psi[87,t,64] <- phiEl[3]*psi4El[9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[87,t,65] <- phiEl[3]*psi4El[9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[87,t,66] <- phiEl[3]*psi4El[9]*(1-p9[6,t])
    for (es in 68:79){
    psi[87,t,es] <- 0}
    psi[87,t,81] <- phiEl[3]*psi4El[10]*p10[6,t]*g10[t, 6]
    psi[87,t,82] <- phiEl[3]*psi4El[10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[87,t,83] <- phiEl[3]*psi4El[10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[87,t,84] <- phiEl[3]*psi4El[10]*(1-p10[6,t])
    psi[87,t,85] <- 0
    psi[87,t,86] <- 0
    psi[87,t,87] <- 0
    psi[87,t,88] <- phiEl[3]*psi4El[11]
    
    for (es in 50:61){
    psi[88,t,es] <- 0}
    psi[88,t,63] <- phiEl[4]*psiAbEl[9]*p9[6,t]*g9[t, 6]
    psi[88,t,64] <- phiEl[4]*psiAbEl[9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[88,t,65] <- phiEl[4]*psiAbEl[9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[88,t,66] <- phiEl[4]*psiAbEl[9]*(1-p9[6,t])
    for (es in 68:79){
    psi[88,t,es] <- 0}
    psi[88,t,81] <- phiEl[4]*psiAbEl[10]*p10[6,t]*g10[t, 6]
    psi[88,t,82] <- phiEl[4]*psiAbEl[10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[88,t,83] <- phiEl[4]*psiAbEl[10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[88,t,84] <- phiEl[4]*psiAbEl[10]*(1-p10[6,t])
    psi[88,t,85] <- 0
    psi[88,t,86] <- 0
    psi[88,t,87] <- 0
    psi[88,t,88] <- phiEl[4]*psiAbEl[11]
    
    
    # Define probabilities of O(t)
    for (j in 1:8){
    for (cc in c(2:4,6)){
    po[((j-1)*6)+cc,t] <- p[j,cc,t] }
    po[1+(6*(j-1)),t] <- 0 #State that can only occur on first capture
    po[5+(6*(j-1)),t] <- 0 #State that can only occur on first capture
    }
    for (jj in c(50:52, 54:56, 58:60, 63:65, 68:70, 72:74, 76:78, 81:83)){
    po[jj,t] <- 1
    }
    po[49,t] <- 0 #State that can only occur on first capture
    po[62,t] <- 0 #State that can only occur on first capture
    po[53,t] <- 0 #Unobservable state
    po[57,t] <- 0 #Unobservable state
    po[61,t] <- 0 #Unobservable state
    po[66,t] <- 0 #Unobservable state
    po[67,t] <- 0 #State that can only occur on first capture
    po[80,t] <- 0 #State that can only occur on first capture
    po[71,t] <- 0 #Unobservable state
    po[75,t] <- 0 #Unobservable state
    po[79,t] <- 0 #Unobservable state
    po[84,t] <- 0 #Unobservable state
    for (j in 85:88){
    po[j,t] <- pEl[j-84]
    }
    # Calculate probability of non-encounter (dq) and reshape the array for the encounter probabilities      
    for (s in 1:ns){
    dp[s,t,s] <- po[s,t]
    dq[s,t,s] <- 1-po[s,t]
    } # s
    
    for (s in 1:(ns-1)){
    for (m in (s+1):ns){
    dp[s,t,m] <- 0
    dq[s,t,m] <- 0
    } # s
    } # m
    for (s in 2:ns){
    for (m in 1:(s-1)){
    dp[s,t,m] <- 0
    dq[s,t,m] <- 0
    } # s
    } # m
    } # t
    
    # Define the multinomial likelihood
    for (t in 1:((n.occasions-1)*ns)){
    marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
    }
    
    # Define the cell probabilities of the multistate m-array   
    # Define matrix U: product of probabilities of state-transition and non-encounter (this is just done because there is no product function for matrix multiplication in JAGS)
    for (t in 1:(n.occasions-2)){
    U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
    for (j in (t+1):(n.occasions-1)){
    U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,t,] %*% dq[,t,]
    }
    }
    U[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ones
    
    # Diagonal
    for (t in 1:(n.occasions-2)){
    pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*% dp[,t,]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% psi[,j,] %*% dp[,j,]
    }
    }
    pr[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- psi[,n.occasions-1,] %*% dp[,n.occasions-1,]
    
    # Below main diagonal
    for (t in 2:(n.occasions-1)){
    for (j in 1:(t-1)){
    pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
    } #j
    } #t
    
    # Last column: probability of non-recapture
    for (t in 1:((n.occasions-1)*ns)){
    pr[t,(n.occasions*ns-(ns-1))] <- 1-sum(pr[t,1:((n.occasions-1)*ns)])
    } #t
    }    ",fill = TRUE)
sink()

# Bundle data
n.unobs <- 8 #Number of unobservable states
n.notobs <- 16 #Number of observable states that dont appear anywhere in capture histories
ns <- length(unique(as.numeric(mat.enc2))) - 1 + n.unobs + n.notobs # calculate the number of states

jags.data <- list(marr = ms.arr, n.occasions = ncol(mat.enc), rel = rowSums(ms.arr), ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), dista=distmat, ind=ind, p.all.years=p.all.years)

inits <- function(){list(a.del=rep(0.5, length.out=3),  phiEl=rep(0.5, length.out=4), 
                         mean.phi = 0.1, sigma.phi = runif(1, 0, 5),
                         mean.p = 0.1, sigma.p = runif(1, 0, 5), A.delp=rnorm(3),
                         t.g9=rnorm(4), a.g.3 = 0.5,a.g.4 = 0.5,a.g.6 = 0.5, t.g10 = 0.5,
                         pEl=c(0.1,0.1,0.1,0.1),
                         B=-0.1,
                         theta2 = rnorm(1), theta3 = rnorm(1), theta4 = rnorm(1), thetaAa = rnorm(1), thetaAb = rnorm(1), repuls=rnorm(10), attract=rnorm(11),
                         mean.psiEl = 0.1)}  

# Parameters monitored
parameters <- c("phi.real", "sigma2.phi.real", "phiEl", "p.real", "sigma2.p.real", "p9.real","g9.real", "p10.real", "g10.real","pEl", "B", "psi1", "psi2","psi3", "psiAa", "psiAb", "psi2El", "psiAbEl", "mean.psiEl")

# MCMC settings
ni <- 1500 #000 
nt <- 1
nb <- 100 #000 
nc <- 3

# Call JAGS from R, using jagsUI syntax
starttime <- Sys.time()
m8 <- jags(jags.data, inits, parameters, "model8.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=FALSE)
endtime <- Sys.time()
options(max.print=999999) #Default is 99999
print(m8, digits = 3) # 3 h for 150 it/chain, 3 chains. B = -0.12
print(m8$summary[1:400,], digits = 3)

#m8first <- m8
#m8second <- m8
quartz()
par(mfrow=c(5,4))
traceplot(m8)

#m8 <- update(m8, n.iter=2000, parallel=TRUE)


#Rerun with different inits
inits.ii <- function(){list(a.del=rep(0.5, length.out=3),  phiEl=rep(0.5, length.out=4), 
                            mean.phi = 0.1, sigma.phi = runif(1, 0, 5),
                            mean.p = 0.1, sigma.p = runif(1, 0, 5), A.delp=rnorm(3),
                            t.g9=rnorm(4), a.g.3 = 0.5,a.g.4 = 0.5,a.g.6 = 0.5, t.g10 = 0.5,
                            pEl=c(0.1,0.1,0.1,0.1),
                            B=rnorm(n=1, mean=-0.1, sd=0.1),
                            theta2 = rnorm(1), theta3 = rnorm(1), theta4 = rnorm(1), thetaAa = rnorm(1), thetaAb = rnorm(1), repuls=rnorm(10), attract=rnorm(11),
                            mean.psiEl = 0.1)}  
# MCMC settings
ni <- 1000 #000 
nt <- 1
nb <- 200 #000 
nc <- 3

starttime <- Sys.time()
m8.ii <- jags(jags.data, inits.ii, parameters, "model8.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
endtime <- Sys.time()
print(m8.ii, digits = 3) #B = -0.12
print(m8.ii$summary[1:400,], digits = 3)

quartz()
par(mfrow=c(5,4))
traceplot(m8.ii)


#Rerun with different inits
inits.iii <- function(){list(a.del=rep(0.5, length.out=3),  phiEl=rep(0.5, length.out=4), 
                             mean.phi = 0.1, sigma.phi = runif(1, 0, 5),
                             mean.p = 0.1, sigma.p = runif(1, 0, 5), A.delp=rnorm(3),
                             t.g9=rnorm(4), a.g.3 = 0.5,a.g.4 = 0.5,a.g.6 = 0.5, t.g10 = 0.5,
                             pEl=c(0.1,0.1,0.1,0.1),
                             B=rnorm(n=1, mean=-0.1, sd=0.2),
                             theta2 = rnorm(1), theta3 = rnorm(1), theta4 = rnorm(1), thetaAa = rnorm(1), thetaAb = rnorm(1), repuls=rnorm(10), attract=rnorm(11),
                             mean.psiEl = 0.1)}  
# MCMC settings
ni <- 500 #000 
nt <- 1
nb <- 200 #000 
nc <- 3

starttime <- Sys.time()
m8.iii <- jags(jags.data, inits.iii, parameters, "model8.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
endtime <- Sys.time()
print(m8.iii, digits = 3) #B = -0.12
print(m8.iii$summary[1:400,], digits = 3)

quartz()
par(mfrow=c(5,4))
traceplot(m8.iii)

#Rerun with different inits
inits.does.not.work <- function(){list(a.del=rep(0.5, length.out=3),  phiEl=rep(0.5, length.out=4), 
                                       mean.phi = 0.1, sigma.phi = runif(1, 0, 5),
                                       mean.p = 0.1, sigma.p = runif(1, 0, 5), A.delp=rnorm(3),
                                       t.g9=rnorm(4), a.g.3 = 0.5,a.g.4 = 0.5,a.g.6 = 0.5, t.g10 = 0.5,
                                       pEl=c(0.1,0.1,0.1,0.1),
                                       B=rnorm(n=1, mean=-0.1, sd=0.5),
                                       theta2 = rnorm(1), theta3 = rnorm(1), theta4 = rnorm(1), thetaAa = rnorm(1), thetaAb = rnorm(1), repuls=rnorm(10), attract=rnorm(11),
                                       mean.psiEl = 0.1)}  

#Rerun with different inits
inits.does.not.work.either <- function(){list(a.del=rep(0.5, length.out=3),  phiEl=rep(0.5, length.out=4), 
                                              mean.phi = 0.1, sigma.phi = runif(1, 0, 5),
                                              mean.p = 0.1, sigma.p = runif(1, 0, 5), A.delp=rnorm(3),
                                              t.g9=rnorm(4), a.g.3 = 0.5,a.g.4 = 0.5,a.g.6 = 0.5, t.g10 = 0.5,
                                              pEl=c(0.1,0.1,0.1,0.1),
                                              B=rnorm(n=1, mean=0, sd=0.1),
                                              theta2 = rnorm(1), theta3 = rnorm(1), theta4 = rnorm(1), thetaAa = rnorm(1), thetaAb = rnorm(1), repuls=rnorm(10), attract=rnorm(11),
                                              mean.psiEl = 0.1)}  

#Rerun with different inits
inits.iv <- function(){list(a.del=rep(0.5, length.out=3),  phiEl=rep(0.5, length.out=4), 
                            mean.phi = 0.1, sigma.phi = runif(1, 0, 5),
                            mean.p = 0.1, sigma.p = runif(1, 0, 5), A.delp=rnorm(3),
                            t.g9=rnorm(4), a.g.3 = 0.5,a.g.4 = 0.5,a.g.6 = 0.5, t.g10 = 0.5,
                            pEl=c(0.1,0.1,0.1,0.1),
                            B=rnorm(n=1, mean=-0.3, sd=0.1),
                            theta2 = rnorm(1), theta3 = rnorm(1), theta4 = rnorm(1), thetaAa = rnorm(1), thetaAb = rnorm(1), repuls=rnorm(10), attract=rnorm(11),
                            mean.psiEl = 0.1)}  

# MCMC settings
ni <- 500 #000 
nt <- 1
nb <- 200 #000 
nc <- 3

starttime <- Sys.time()
m8.iv <- jags(jags.data, inits.iv, parameters, "model8.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
endtime <- Sys.time()
print(m8.iv, digits = 3) #B = -0.12
print(m8.iv$summary[1:400,], digits = 3)

quartz()
par(mfrow=c(5,4))
traceplot(m8.iv)









########################
########################
#Model 9
########################
#Removing the Elsewhere states (and of course transtions to and from them, and encounter probabilities)

#Creating an index matrix, to use in script
in. <- matrix(data=rep(1:9,10), ncol=9, nrow=10, byrow=T)
for (k in 1:9){in.[k,k] <- 10}

#Study characteristics
n.occasions <- 10 #Using the data from the last ten years only 
#(thus some birds will appear as they were first marked as 2/3/4-years-olds etc)
n.states <- 85 #Removing the 4 Elsewhere states
n.obs <- 85 

#Removing the first years of data
enc0817 <- enc[, c(1:5, 19:31)] 

#Combining the numbers for adults as chicks as information about ringing stage is incorporated in the state
enc0817$No. <- enc0817$Ch + enc0817$Ad
#Selecting the relevant columns and changing the order
enc0817. <- enc0817[, c(18, 2:15, 19)] 

#Creating one encounter history for each individual
mat.enc <- matrix(NA, ncol = n.occasions, nrow = sum(enc0817.[,ncol(enc0817.)]))
i <- 1
for (l in 1:enc0817.[i,ncol(enc0817.)]){
  for (k in 1:n.occasions){
    mat.enc[l,k] <- enc0817.[i, (ncol(enc0817.)-n.occasions+k-1)]  } }

for (i in 2:nrow(enc0817.)){
  j <- 1 + sum(enc0817.[1:(i-1),ncol(enc0817.)])
  for (l in j:(sum(enc0817.[1:(i-1),ncol(enc0817.)])+enc0817.[i,ncol(enc0817.)])){
    for (k in 1:n.occasions){
      mat.enc[l,k] <- enc0817.[i, (ncol(enc0817.)-n.occasions+k-1)]    } }}

#Removing the rows with no events (only 0)
mat.enc2 <- mat.enc[which(rowSums(mat.enc)>0),]

#Checking the encounter histories
table(mat.enc2) #Need to check if there are any events that didn't happen, as some probabilities are low

#Recoding the Elsewhere states to 0
elsew <- which(mat.enc2==88|mat.enc2==85|mat.enc2==86)
mat.enc3 <- mat.enc2
mat.enc3[elsew] <- 0

#Checking rows with no events (only 0) 
mat.enc3[which(rowSums(mat.enc3)==0),]
mat.enc2[which(rowSums(mat.enc3)==0),]
#Removing the rows with no events (only 0) - individuals marked prior to first year analyzed, and only encountered Elsewhere, hence with only 0 in encounter history
mat.enc3 <- mat.enc3[which(rowSums(mat.enc3)>0),]
# Compute vector with occasion of first capture
f <- apply(mat.enc3, 1, get.first)

#Checking the encounter histories
table(mat.enc3) 
length(table(mat.enc3))
#Reformating the data using the m-array function. 
ms.arr <- marray(mat.enc3, unobs=23) #Unobs needs to be no. of unobservable states + states not observed (despite observable)

p.all.years <- c(1,3,8) #Sites 1,3,8 should have p estimated for all years of the study

# Analysis of the model
# Analysis of the model
# Specify model in BUGS language
sink("model9.jags") 
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi[fro, stage, t]: survival probability at age/stage, having been at site fro, in period t
    # psi[fro, till]: movement probability from site fro to site till, only 10 of these as stage at next time-step is known and thereby constraints options to the 10 sites
    # p[...]: encounter probability. 
    # g: proportion of birds encountered (observed, captured or observed & captured) that where captured
    # -------------------------------------------------
    # States (S):
    # 1 alive at a1 s1
    # 2 alive at a2 s1
    # 3 alive at a3 s1
    # 4 alive at a4 s1
    # 5 alive at Ad_recent s1
    # 6 alive at Ad_previous s1
    # 7 alive at a1 s2
    # ...
    # 48 alive at Ad_previous s8
    # 49 alive & captured at a1 s9
    # 50 alive & captured at a2 s9
    # 51 alive & observed at a2 s9
    # 52 alive & observed & captured at a2 s9
    # 53 alive but not seen at a2 s9
    # 54 alive & captured at a3 s9
    # ...
    # 84 alive but not seen at Ad_previous s10
    # 85 dead
    #
    # Observations (O):
    # 1 seen at a1 s1
    # 2 seen at a2 s1
    # 3 seen at a3 s1
    # 4 seen at a4 s1
    # 5 seen at Ad s1
    # 6 seen at a1 s2
    # ...
    # 85/0 not seen
    # Unobservable states: 53,57,61,66, 71,75,79,84
    
    # -------------------------------------------------
    
    # Priors and constraints
    # Survival: random effects of time and island, additive effect of age
    for (fro in 1:10){
    for (stage in 1:6){ 
    for (t in 1:(n.occasions-1)){
    logit(phi[fro, stage, t]) <- mu.phi + eps.phi[fro, t] + del[stage] 
    }}}
    for (fro in 1:10){
    for (t in 1:(n.occasions-1)){
    eps.phi[fro, t] ~ dnorm(0, tau.phi)}}
    mean.phi ~ dunif(0, 1)                    # Prior for mean survival
    mu.phi <- log(mean.phi / (1-mean.phi))    # Logit transformation
    sigma.phi ~ dunif(0, 5)                       # Prior for standard deviation
    tau.phi <- pow(sigma.phi, -2)
    sigma2.phi <- pow(sigma.phi, 2)
    sigma2.phi.real <- sigma2.phi * pow(mean.phi, 2) * pow((1-mean.phi), 2) # Temporal variance on real scale
    
    del[1] <- 0                             # Priors for age/stage-effect
    for (stag in 2:4){ 
    del[stag] <- a.del[1]}
    a.del[1] ~ dnorm(0, 0.01)
    for (stag in 5:6){ 
    del[stag] <- a.del[stag-3]
    a.del[stag-3] ~ dnorm(0, 0.01)}
    for (fron in 1:10){  #Back-transform to the probability scale
    for (stage in 1:6){ 
    for (t in 1:(n.occasions-1)){
    phi.real[fron, stage, t] <- 1/(1+exp(-(mu.phi + eps.phi[fron,t] + del[stage])))    }}}
    
    # Encounters at site 1-8: random effect of time and island, additive effect of age
    #Sites 1,3,8: 
    for (isl in p.all.years){ 
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[isl, stage, t]) <- mu.p + eps.p[isl, t] + delp[stage]     }}}
    for (isl in p.all.years){
    for (t in 1:(n.occasions-1)){
    eps.p[isl, t] ~ dnorm(0, tau.p)}}
    for (isln in p.all.years){  #Back-transform to the probability scale
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    p.real[isln, stage, t] <- 1/(1+exp(-(mu.p + eps.p[isln,t] + delp[stage])))    }}}
    #Site 2 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[2, stage, t]) <- mu.p + eps.p[2, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[2, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[2, stage, t] <- 1/(1+exp(-(mu.p + eps.p[2,t] + delp[stage])))    }}
    #Site 4 (should have p estimated for 2009-2011,2013-2017):
    for (stage in c(2:4,6)){ 
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    logit(p[4, stage, t]) <- mu.p + eps.p[4, t] + delp[stage]   }}
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    eps.p[4, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    p.real[4, stage, t] <- 1/(1+exp(-(mu.p + eps.p[4,t] + delp[stage])))    }}
    for (stage in c(2:4,6)){ 
    for (t in (n.occasions-6)){
    p[4, stage, t] <- 0   }}
    #Site 5 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[5, stage, t]) <- mu.p + eps.p[5, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[5, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[5, stage, t] <- 1/(1+exp(-(mu.p + eps.p[5,t] + delp[stage])))    }}
    #Site 6 (should have p estimated for 2011-2017):
    for (stage in c(2:4,6)){ 
    for (t in (n.occasions-7):(n.occasions-1)){
    logit(p[6, stage, t]) <- mu.p + eps.p[6, t] + delp[stage]   }}
    for (t in (n.occasions-7):(n.occasions-1)){
    eps.p[6, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in (n.occasions-7):(n.occasions-1)){
    p.real[6, stage, t] <- 1/(1+exp(-(mu.p + eps.p[6,t] + delp[stage])))    }}
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-8)){
    p[6, stage, t] <- 0   }}
    #Site 7 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[7, stage, t]) <- mu.p + eps.p[7, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[7, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[7, stage, t] <- 1/(1+exp(-(mu.p + eps.p[7,t] + delp[stage])))    }}
    
    mean.p ~ dunif(0, 1)                    # Prior for mean re-encounter
    mu.p <- log(mean.p / (1-mean.p))        # Logit transformation
    sigma.p ~ dunif(0, 5)                   # Prior for standard deviation
    tau.p <- pow(sigma.p, -2)
    sigma2.p <- pow(sigma.p, 2)
    sigma2.p.real <- sigma2.p * pow(mean.p, 2) * pow((1-mean.p), 2) # Temporal variance on real scale
    #delp[1] <- 0 #Re-encounters do not happen for age 1 [Shouldn't be necessary to define when using stage c(2:4,6) everywhere for p]
    #delp[5] <- 0 #Re-encounters do not happen for full-grown, just marked [Shouldn't be necessary to define when using stage c(2:4,6) everywhere for p]
    delp[2] <- 0 #Baseline age-class, which mu.p applies to
    for (stag in c(3,4)){ # Priors for age/stage-effect
    delp[stag] <- A.delp[stag-2]
    A.delp[stag-2] ~ dnorm(0, 0.01)}
    delp[6] <- A.delp[3]
    A.delp[3] ~ dnorm(0, 0.01)
    
    # Encounters at site 9-10:
    #Site 9 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p9[stage, t]) <- mu.p + eps.p9[t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p9[t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p9.real[stage, t] <- 1/(1+exp(-(mu.p + eps.p9[t] + delp[stage])))    }}
    #Site 10 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p10[stage, t]) <- mu.p + eps.p10[t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p10[t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p10.real[stage, t] <- 1/(1+exp(-(mu.p + eps.p10[t] + delp[stage])))    }}
    
    for (fd in 1:(n.occasions-5)){
    for (st in c(2:4,6)){
    g9[fd, st] <- 1  }}
    for (fd in (n.occasions-4):(n.occasions-1)){
    logit(g9[fd, 2]) <- t.g9[fd-n.occasions+5]
    logit(g9[fd, 3]) <- t.g9[fd-n.occasions+5] + a.g.3
    logit(g9[fd, 4]) <- t.g9[fd-n.occasions+5] + a.g.4
    logit(g9[fd, 6]) <- t.g9[fd-n.occasions+5] + a.g.6
    t.g9[fd-n.occasions+5] ~ dnorm(0, 0.01)
    g9.real[fd,2] <- 1/(1+exp(-(t.g9[fd-n.occasions+5]))) 
    g9.real[fd,3] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.3)))
    g9.real[fd,4] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.4)))
    g9.real[fd,6] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.6)))}
    for (fd in 1:(n.occasions-4)){ 
    for (st in c(2:4,6)){
    g10[fd, st] <- 1  }}
    logit(g10[n.occasions-3, 2]) <- t.g10
    logit(g10[n.occasions-3, 3]) <- t.g10 + a.g.3
    logit(g10[n.occasions-3, 4]) <- t.g10 + a.g.4
    logit(g10[n.occasions-3, 6]) <- t.g10 + a.g.6
    t.g10 ~ dnorm(0, 0.01)
    g10.real[n.occasions-3,2] <- 1/(1+exp(-(t.g10))) 
    g10.real[n.occasions-3,3] <- 1/(1+exp(-(t.g10 + a.g.3)))
    g10.real[n.occasions-3,4] <- 1/(1+exp(-(t.g10 + a.g.4)))
    g10.real[n.occasions-3,6] <- 1/(1+exp(-(t.g10 + a.g.6)))
    for (fd in (n.occasions-2):(n.occasions-1)){ 
    for (st in c(2:4,6)){
    g10[fd, st] <- 1  }}
    a.g.3 ~ dnorm(0, 0.01)
    a.g.4 ~ dnorm(0, 0.01)
    a.g.6 ~ dnorm(0, 0.01)
    
    # Transitions: multinomial logit
    # Normal priors on logit of all but one transition probs
    for (fro in 1:10){
    for (till in 1:(fro-1)){
    lpsi1[fro, till] <- repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi2[fro, till] <- theta2 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi3[fro, till] <- theta3 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi4[fro, till] <- theta4 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAa[fro, till] <- thetaAa + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAb[fro, till] <- thetaAb + repuls[fro] + attract[till] + B*dista[fro,till]}
    for (till in (fro+1):10){
    lpsi1[fro, till] <- repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi2[fro, till] <- theta2 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi3[fro, till] <- theta3 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi4[fro, till] <- theta4 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAa[fro, till] <- thetaAa + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAb[fro, till] <- thetaAb + repuls[fro] + attract[till] + B*dista[fro,till]} }
    
    for (fro in 1:10){
    for (till in 1:(fro-1)){
    psi1[fro, till] <- exp(lpsi1[fro, till]) / (1 + exp(lpsi1[fro,ind[fro,1]]) + exp(lpsi1[fro,ind[fro,2]]) + exp(lpsi1[fro,ind[fro,3]]) + exp(lpsi1[fro,ind[fro,4]]) + exp(lpsi1[fro,ind[fro,5]]) + exp(lpsi1[fro,ind[fro,6]]) + exp(lpsi1[fro,ind[fro,7]]) + exp(lpsi1[fro,ind[fro,8]]) + exp(lpsi1[fro,ind[fro,9]]))
    psi2[fro, till] <- exp(lpsi2[fro, till]) / (1 + exp(lpsi2[fro,ind[fro,1]]) + exp(lpsi2[fro,ind[fro,2]]) + exp(lpsi2[fro,ind[fro,3]]) + exp(lpsi2[fro,ind[fro,4]]) + exp(lpsi2[fro,ind[fro,5]]) + exp(lpsi2[fro,ind[fro,6]]) + exp(lpsi2[fro,ind[fro,7]]) + exp(lpsi2[fro,ind[fro,8]]) + exp(lpsi2[fro,ind[fro,9]]))
    psi3[fro, till] <- exp(lpsi3[fro, till]) / (1 + exp(lpsi3[fro,ind[fro,1]]) + exp(lpsi3[fro,ind[fro,2]]) + exp(lpsi3[fro,ind[fro,3]]) + exp(lpsi3[fro,ind[fro,4]]) + exp(lpsi3[fro,ind[fro,5]]) + exp(lpsi3[fro,ind[fro,6]]) + exp(lpsi3[fro,ind[fro,7]]) + exp(lpsi3[fro,ind[fro,8]]) + exp(lpsi3[fro,ind[fro,9]]))
    psi4[fro, till] <- exp(lpsi4[fro, till]) / (1 + exp(lpsi4[fro,ind[fro,1]]) + exp(lpsi4[fro,ind[fro,2]]) + exp(lpsi4[fro,ind[fro,3]]) + exp(lpsi4[fro,ind[fro,4]]) + exp(lpsi4[fro,ind[fro,5]]) + exp(lpsi4[fro,ind[fro,6]]) + exp(lpsi4[fro,ind[fro,7]]) + exp(lpsi4[fro,ind[fro,8]]) + exp(lpsi4[fro,ind[fro,9]]))
    psiAa[fro, till] <- exp(lpsiAa[fro, till]) / (1 + exp(lpsiAa[fro,ind[fro,1]]) + exp(lpsiAa[fro,ind[fro,2]]) + exp(lpsiAa[fro,ind[fro,3]]) + exp(lpsiAa[fro,ind[fro,4]]) + exp(lpsiAa[fro,ind[fro,5]]) + exp(lpsiAa[fro,ind[fro,6]]) + exp(lpsiAa[fro,ind[fro,7]]) + exp(lpsiAa[fro,ind[fro,8]]) + exp(lpsiAa[fro,ind[fro,9]]))
    psiAb[fro, till] <- exp(lpsiAb[fro, till]) / (1 + exp(lpsiAb[fro,ind[fro,1]]) + exp(lpsiAb[fro,ind[fro,2]]) + exp(lpsiAb[fro,ind[fro,3]]) + exp(lpsiAb[fro,ind[fro,4]]) + exp(lpsiAb[fro,ind[fro,5]]) + exp(lpsiAb[fro,ind[fro,6]]) + exp(lpsiAb[fro,ind[fro,7]]) + exp(lpsiAb[fro,ind[fro,8]]) + exp(lpsiAb[fro,ind[fro,9]]))  }
    for (till in (fro+1):10){
    psi1[fro, till] <- exp(lpsi1[fro, till]) / (1 + exp(lpsi1[fro,ind[fro,1]]) + exp(lpsi1[fro,ind[fro,2]]) + exp(lpsi1[fro,ind[fro,3]]) + exp(lpsi1[fro,ind[fro,4]]) + exp(lpsi1[fro,ind[fro,5]]) + exp(lpsi1[fro,ind[fro,6]]) + exp(lpsi1[fro,ind[fro,7]]) + exp(lpsi1[fro,ind[fro,8]]) + exp(lpsi1[fro,ind[fro,9]]))
    psi2[fro, till] <- exp(lpsi2[fro, till]) / (1 + exp(lpsi2[fro,ind[fro,1]]) + exp(lpsi2[fro,ind[fro,2]]) + exp(lpsi2[fro,ind[fro,3]]) + exp(lpsi2[fro,ind[fro,4]]) + exp(lpsi2[fro,ind[fro,5]]) + exp(lpsi2[fro,ind[fro,6]]) + exp(lpsi2[fro,ind[fro,7]]) + exp(lpsi2[fro,ind[fro,8]]) + exp(lpsi2[fro,ind[fro,9]]))
    psi3[fro, till] <- exp(lpsi3[fro, till]) / (1 + exp(lpsi3[fro,ind[fro,1]]) + exp(lpsi3[fro,ind[fro,2]]) + exp(lpsi3[fro,ind[fro,3]]) + exp(lpsi3[fro,ind[fro,4]]) + exp(lpsi3[fro,ind[fro,5]]) + exp(lpsi3[fro,ind[fro,6]]) + exp(lpsi3[fro,ind[fro,7]]) + exp(lpsi3[fro,ind[fro,8]]) + exp(lpsi3[fro,ind[fro,9]]))
    psi4[fro, till] <- exp(lpsi4[fro, till]) / (1 + exp(lpsi4[fro,ind[fro,1]]) + exp(lpsi4[fro,ind[fro,2]]) + exp(lpsi4[fro,ind[fro,3]]) + exp(lpsi4[fro,ind[fro,4]]) + exp(lpsi4[fro,ind[fro,5]]) + exp(lpsi4[fro,ind[fro,6]]) + exp(lpsi4[fro,ind[fro,7]]) + exp(lpsi4[fro,ind[fro,8]]) + exp(lpsi4[fro,ind[fro,9]]))
    psiAa[fro, till] <- exp(lpsiAa[fro, till]) / (1 + exp(lpsiAa[fro,ind[fro,1]]) + exp(lpsiAa[fro,ind[fro,2]]) + exp(lpsiAa[fro,ind[fro,3]]) + exp(lpsiAa[fro,ind[fro,4]]) + exp(lpsiAa[fro,ind[fro,5]]) + exp(lpsiAa[fro,ind[fro,6]]) + exp(lpsiAa[fro,ind[fro,7]]) + exp(lpsiAa[fro,ind[fro,8]]) + exp(lpsiAa[fro,ind[fro,9]]))
    psiAb[fro, till] <- exp(lpsiAb[fro, till]) / (1 + exp(lpsiAb[fro,ind[fro,1]]) + exp(lpsiAb[fro,ind[fro,2]]) + exp(lpsiAb[fro,ind[fro,3]]) + exp(lpsiAb[fro,ind[fro,4]]) + exp(lpsiAb[fro,ind[fro,5]]) + exp(lpsiAb[fro,ind[fro,6]]) + exp(lpsiAb[fro,ind[fro,7]]) + exp(lpsiAb[fro,ind[fro,8]]) + exp(lpsiAb[fro,ind[fro,9]]))  }
    #The case below is when fro=till, but rewritten to match JAGS syntax
    psi1[fro,fro] <- 1- sum(psi1[fro,ind[fro,]]) 
    psi2[fro,fro] <- 1- sum(psi2[fro,ind[fro,]]) 
    psi3[fro,fro] <- 1- sum(psi3[fro,ind[fro,]]) 
    psi4[fro,fro] <- 1- sum(psi4[fro,ind[fro,]]) 
    psiAa[fro,fro] <- 1- sum(psiAa[fro,ind[fro,]])   
    psiAb[fro,fro] <- 1- sum(psiAb[fro,ind[fro,]])   
    } 
    
    #Priors for transitions
    theta2 ~ dnorm(0, 0.01)               # Prior for age-effect in psi
    theta3 ~ dnorm(0, 0.01)  
    theta4 ~ dnorm(0, 0.01)  
    thetaAa ~ dnorm(0, 0.01)  
    thetaAb ~ dnorm(0, 0.01) 
    #mean.psi ~ dunif(0, 1)                # Prior for mean movement prob - not used as no transition makes sense to set as baseline
    #mu <- log(mean.psi / (1-mean.psi))    # Logit transformation
    for (w in 1:10){
    repuls[w] ~ dnorm(0, 0.01) 
    attract[w] ~ dnorm(0, 0.01) 
    }
    B ~ dnorm(0, 0.001)I(-10,10)
    
    # Define state-transition and observation matrices 	
    # Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
    for (d in 1:84){
    for (m in 1:8){
    psi[d,t,1+(6*(m-1))] <- 0
    psi[d,t,5+(6*(m-1))] <- 0
    } #m
    psi[d,t,49] <- 0
    psi[d,t,62] <- 0
    psi[d,t,67] <- 0
    psi[d,t,80] <- 0
    } #d
    
    for (sf in 1:8){  
    for (s in 1:8){
    psi[(6*(sf-1))+1,t,(s*6-4)] <- phi[sf,1,t]*psi1[sf,s]
    psi[(6*(sf-1))+1,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+1,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+1,t,(s*6)] <- 0
    psi[(6*(sf-1))+2,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+2,t,(s*6-3)] <- phi[sf,2,t]*psi2[sf,s]
    psi[(6*(sf-1))+2,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+2,t,(s*6)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-2)] <- phi[sf,3,t]*psi3[sf,s]
    psi[(6*(sf-1))+3,t,(s*6)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+4,t,(s*6)] <- phi[sf,4,t]*psi4[sf,s]
    psi[(6*(sf-1))+5,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+5,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+5,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+5,t,(s*6)] <- phi[sf,5,t]*psiAa[sf,s]
    psi[(6*(sf-1))+6,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+6,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+6,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+6,t,(s*6)] <- phi[sf,6,t]*psiAb[sf,s]
    } #s
    
    psi[(6*(sf-1))+1,t,50] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*g9[t, 2]
    psi[(6*(sf-1))+1,t,51] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[(6*(sf-1))+1,t,52] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[(6*(sf-1))+1,t,53] <- phi[sf,1,t]*psi1[sf,9]*(1-p9[2,t])
    for (ex in c(54:61,63:66)){
    psi[(6*(sf-1))+1,t,ex] <- 0  }
    psi[(6*(sf-1))+1,t,68] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*g10[t, 2]
    psi[(6*(sf-1))+1,t,69] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[(6*(sf-1))+1,t,70] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[(6*(sf-1))+1,t,71] <- phi[sf,1,t]*psi1[sf,10]*(1-p10[2,t])
    for (ex in c(72:79,81:84)){
    psi[(6*(sf-1))+1,t,ex] <- 0  }
    
    for (ew in 50:53){
    psi[(6*(sf-1))+2,t,ew] <- 0}
    psi[(6*(sf-1))+2,t,54] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*g9[t, 3]
    psi[(6*(sf-1))+2,t,55] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[(6*(sf-1))+2,t,56] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[(6*(sf-1))+2,t,57] <- phi[sf,2,t]*psi2[sf,9]*(1-p9[3,t])
    for (er in c(58:61,63:66)){
    psi[(6*(sf-1))+2,t,er] <- 0  }
    for (ew in 68:71){
    psi[(6*(sf-1))+2,t,ew] <- 0}
    psi[(6*(sf-1))+2,t,72] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*g10[t, 3]
    psi[(6*(sf-1))+2,t,73] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[(6*(sf-1))+2,t,74] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[(6*(sf-1))+2,t,75] <- phi[sf,2,t]*psi2[sf,10]*(1-p10[3,t])
    for (er in c(76:79,81:84)){
    psi[(6*(sf-1))+2,t,er] <- 0  }
    
    for (ed in 50:57){
    psi[(6*(sf-1))+3,t,ed] <- 0}
    psi[(6*(sf-1))+3,t,58] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*g9[t, 4]
    psi[(6*(sf-1))+3,t,59] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[(6*(sf-1))+3,t,60] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[(6*(sf-1))+3,t,61] <- phi[sf,3,t]*psi3[sf,9]*(1-p9[4,t])
    for (ef in 63:66){
    psi[(6*(sf-1))+3,t,ef] <- 0 }
    for (ed in 68:75){
    psi[(6*(sf-1))+3,t,ed] <- 0}
    psi[(6*(sf-1))+3,t,76] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*g10[t, 4]
    psi[(6*(sf-1))+3,t,77] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[(6*(sf-1))+3,t,78] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[(6*(sf-1))+3,t,79] <- phi[sf,3,t]*psi3[sf,10]*(1-p10[4,t])
    for (ef in 81:84){
    psi[(6*(sf-1))+3,t,ef] <- 0 }
    
    for (es in 50:61){
    psi[(6*(sf-1))+4,t,es] <- 0}
    psi[(6*(sf-1))+4,t,63] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+4,t,64] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+4,t,65] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+4,t,66] <- phi[sf,4,t]*psi4[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+4,t,es] <- 0}
    psi[(6*(sf-1))+4,t,81] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+4,t,82] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+4,t,83] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+4,t,84] <- phi[sf,4,t]*psi4[sf,10]*(1-p10[6,t])
    
    for (es in 50:61){
    psi[(6*(sf-1))+5,t,es] <- 0}
    psi[(6*(sf-1))+5,t,63] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+5,t,64] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+5,t,65] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+5,t,66] <- phi[sf,5,t]*psiAa[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+5,t,es] <- 0}
    psi[(6*(sf-1))+5,t,81] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+5,t,82] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+5,t,83] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+5,t,84] <- phi[sf,5,t]*psiAa[sf,10]*(1-p10[6,t])
    
    for (es in 50:61){
    psi[(6*(sf-1))+6,t,es] <- 0}
    psi[(6*(sf-1))+6,t,63] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+6,t,64] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+6,t,65] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+6,t,66] <- phi[sf,6,t]*psiAb[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+6,t,es] <- 0}
    psi[(6*(sf-1))+6,t,81] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+6,t,82] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+6,t,83] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+6,t,84] <- phi[sf,6,t]*psiAb[sf,10]*(1-p10[6,t])
    } #sf
    
    
    for (s in 1:8){
    psi[49,t,(s*6-4)] <- phi[9,1,t]*psi1[9,s]
    psi[49,t,(s*6-3)] <- 0
    psi[49,t,(s*6-2)] <- 0
    psi[49,t,(s*6)] <- 0
    for (i in 1:4){
    psi[49+i,t,(s*6-4)] <- 0
    psi[49+i,t,(s*6-3)] <- phi[9,2,t]*psi2[9,s]
    psi[49+i,t,(s*6-2)] <- 0
    psi[49+i,t,(s*6)] <- 0
    psi[53+i,t,(s*6-4)] <- 0
    psi[53+i,t,(s*6-3)] <- 0
    psi[53+i,t,(s*6-2)] <- phi[9,3,t]*psi3[9,s]
    psi[53+i,t,(s*6)] <- 0
    psi[57+i,t,(s*6-4)] <- 0
    psi[57+i,t,(s*6-3)] <- 0
    psi[57+i,t,(s*6-2)] <- 0
    psi[57+i,t,(s*6)] <- phi[9,4,t]*psi4[9,s] 
    psi[62+i,t,(s*6-4)] <- 0
    psi[62+i,t,(s*6-3)] <- 0
    psi[62+i,t,(s*6-2)] <- 0
    psi[62+i,t,(s*6)] <- phi[9,6,t]*psiAb[9,s]
    } #i
    psi[62,t,(s*6-4)] <- 0
    psi[62,t,(s*6-3)] <- 0
    psi[62,t,(s*6-2)] <- 0
    psi[62,t,(s*6)] <- phi[9,5,t]*psiAa[9,s]
    } #s
    
    for (s in 1:8){
    psi[67,t,(s*6-4)] <- phi[10,1,t]*psi1[10,s]
    psi[67,t,(s*6-3)] <- 0
    psi[67,t,(s*6-2)] <- 0
    psi[67,t,(s*6)] <- 0
    for (i in 1:4){
    psi[67+i,t,(s*6-4)] <- 0
    psi[67+i,t,(s*6-3)] <- phi[10,2,t]*psi2[10,s]
    psi[67+i,t,(s*6-2)] <- 0
    psi[67+i,t,(s*6)] <- 0
    psi[71+i,t,(s*6-4)] <- 0
    psi[71+i,t,(s*6-3)] <- 0
    psi[71+i,t,(s*6-2)] <- phi[10,3,t]*psi3[10,s]
    psi[71+i,t,(s*6)] <- 0
    psi[75+i,t,(s*6-4)] <- 0
    psi[75+i,t,(s*6-3)] <- 0
    psi[75+i,t,(s*6-2)] <- 0
    psi[75+i,t,(s*6)] <- phi[10,4,t]*psi4[10,s] 
    psi[80+i,t,(s*6-4)] <- 0
    psi[80+i,t,(s*6-3)] <- 0
    psi[80+i,t,(s*6-2)] <- 0
    psi[80+i,t,(s*6)] <- phi[10,6,t]*psiAb[10,s]
    } #i
    psi[80,t,(s*6-4)] <- 0
    psi[80,t,(s*6-3)] <- 0
    psi[80,t,(s*6-2)] <- 0
    psi[80,t,(s*6)] <- phi[10,5,t]*psiAa[10,s]
    } #s
    
    psi[49,t,50] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*g9[t, 2]
    psi[49,t,51] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[49,t,52] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[49,t,53] <- phi[9,1,t]*psi1[9,9]*(1-p9[2,t])
    for (ex in c(54:61,63:66)){
    psi[49,t,ex] <- 0  }
    
    for (i in 1:4){
    for (ew in 50:53){
    psi[49+i,t,ew] <- 0  }
    psi[49+i,t,54] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*g9[t, 3]
    psi[49+i,t,55] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[49+i,t,56] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[49+i,t,57] <- phi[9,2,t]*psi2[9,9]*(1-p9[3,t])
    for (er in c(58:61,63:66)){
    psi[49+i,t,er] <- 0  }
    for (ed in 50:57){
    psi[53+i,t,ed] <- 0}
    psi[53+i,t,58] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*g9[t, 4]
    psi[53+i,t,59] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[53+i,t,60] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[53+i,t,61] <- phi[9,3,t]*psi3[9,9]*(1-p9[4,t])
    for (ef in 63:66){
    psi[53+i,t,ef] <- 0 }
    for (es in 50:61){
    psi[57+i,t,es] <- 0}
    psi[57+i,t,63] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*g9[t, 6]
    psi[57+i,t,64] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[57+i,t,65] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[57+i,t,66] <- phi[9,4,t]*psi4[9,9]*(1-p9[6,t])
    for (es in 50:61){
    psi[62+i,t,es] <- 0}
    psi[62+i,t,63] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*g9[t, 6]
    psi[62+i,t,64] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62+i,t,65] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62+i,t,66] <- phi[9,6,t]*psiAb[9,9]*(1-p9[6,t])
    } #i
    
    for (es in 50:61){
    psi[62,t,es] <- 0}
    psi[62,t,63] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*g9[t, 6]
    psi[62,t,64] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62,t,65] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62,t,66] <- phi[9,5,t]*psiAa[9,9]*(1-p9[6,t])
    
    
    psi[49,t,50+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*g10[t, 2]
    psi[49,t,51+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[49,t,52+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[49,t,53+18] <- phi[9,1,t]*psi1[9,10]*(1-p10[2,t])
    for (ex in c(54:61,63:66)){
    psi[49,t,ex+18] <- 0  }
    
    for (i in 1:4){
    for (ew in 50:53){
    psi[49+i,t,ew+18] <- 0  }
    psi[49+i,t,54+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*g10[t, 3]
    psi[49+i,t,55+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[49+i,t,56+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[49+i,t,57+18] <- phi[9,2,t]*psi2[9,10]*(1-p10[3,t])
    for (er in c(58:61,63:66)){
    psi[49+i,t,er+18] <- 0  }
    for (ed in 50:57){
    psi[53+i,t,ed+18] <- 0}
    psi[53+i,t,58+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*g10[t, 4]
    psi[53+i,t,59+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[53+i,t,60+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[53+i,t,61+18] <- phi[9,3,t]*psi3[9,10]*(1-p10[4,t])
    for (ef in 63:66){
    psi[53+i,t,ef+18] <- 0 }
    for (es in 50:61){
    psi[57+i,t,es+18] <- 0}
    psi[57+i,t,63+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*g10[t, 6]
    psi[57+i,t,64+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[57+i,t,65+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[57+i,t,66+18] <- phi[9,4,t]*psi4[9,10]*(1-p10[6,t])
    for (es in 50:61){
    psi[62+i,t,es+18] <- 0}
    psi[62+i,t,63+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*g10[t, 6]
    psi[62+i,t,64+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62+i,t,65+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62+i,t,66+18] <- phi[9,6,t]*psiAb[9,10]*(1-p10[6,t])
    } #i
    
    for (es in 50:61){
    psi[62,t,es+18] <- 0}
    psi[62,t,63+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*g10[t, 6]
    psi[62,t,64+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62,t,65+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62,t,66+18] <- phi[9,5,t]*psiAa[9,10]*(1-p10[6,t])
    
    psi[67,t,50] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*g9[t, 2]
    psi[67,t,51] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[67,t,52] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[67,t,53] <- phi[10,1,t]*psi1[10,9]*(1-p9[2,t])
    for (ex in c(72:79,81:84)){
    psi[67,t,ex-18] <- 0  }
    
    for (i in 1:4){
    for (ew in 68:71){
    psi[67+i,t,ew-18] <- 0  }
    psi[67+i,t,72-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*g9[t, 3]
    psi[67+i,t,73-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[67+i,t,74-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[67+i,t,75-18] <- phi[10,2,t]*psi2[10,9]*(1-p9[3,t])
    for (er in c(76:79,81:84)){
    psi[67+i,t,er-18] <- 0  }
    for (ed in 68:75){
    psi[71+i,t,ed-18] <- 0}
    psi[71+i,t,76-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*g9[t, 4]
    psi[71+i,t,77-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[71+i,t,78-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[71+i,t,79-18] <- phi[10,3,t]*psi3[10,9]*(1-p9[4,t])
    for (ef in 81:84){
    psi[71+i,t,ef-18] <- 0 }
    for (es in 68:79){
    psi[75+i,t,es-18] <- 0}
    psi[75+i,t,81-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*g9[t, 6]
    psi[75+i,t,82-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[75+i,t,83-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[75+i,t,84-18] <- phi[10,4,t]*psi4[10,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[80+i,t,es-18] <- 0}
    psi[80+i,t,81-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*g9[t, 6]
    psi[80+i,t,82-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80+i,t,83-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80+i,t,84-18] <- phi[10,6,t]*psiAb[10,9]*(1-p9[6,t])
    } #i
    
    for (es in 68:79){
    psi[80,t,es-18] <- 0}
    psi[80,t,81-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*g9[t, 6]
    psi[80,t,82-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80,t,83-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80,t,84-18] <- phi[10,5,t]*psiAa[10,9]*(1-p9[6,t])
    
    
    psi[67,t,68] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*g10[t, 2]
    psi[67,t,69] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[67,t,70] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[67,t,71] <- phi[10,1,t]*psi1[10,10]*(1-p10[2,t])
    for (ex in c(72:79,81:84)){
    psi[67,t,ex] <- 0  }
    
    for (i in 1:4){
    for (ew in 68:71){
    psi[67+i,t,ew] <- 0  }
    psi[67+i,t,72] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*g10[t, 3]
    psi[67+i,t,73] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[67+i,t,74] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[67+i,t,75] <- phi[10,2,t]*psi2[10,10]*(1-p10[3,t])
    for (er in c(76:79,81:84)){
    psi[67+i,t,er] <- 0  }
    for (ed in 68:75){
    psi[71+i,t,ed] <- 0}
    psi[71+i,t,76] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*g10[t, 4]
    psi[71+i,t,77] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[71+i,t,78] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[71+i,t,79] <- phi[10,3,t]*psi3[10,10]*(1-p10[4,t])
    for (ef in 81:84){
    psi[71+i,t,ef] <- 0 }
    for (es in 68:79){
    psi[75+i,t,es] <- 0}
    psi[75+i,t,81] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*g10[t, 6]
    psi[75+i,t,82] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[75+i,t,83] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[75+i,t,84] <- phi[10,4,t]*psi4[10,10]*(1-p10[6,t])
    for (es in 68:79){
    psi[80+i,t,es] <- 0}
    psi[80+i,t,81] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*g10[t, 6]
    psi[80+i,t,82] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80+i,t,83] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80+i,t,84] <- phi[10,6,t]*psiAb[10,10]*(1-p10[6,t])
    } #i
    for (es in 68:79){
    psi[80,t,es] <- 0}
    psi[80,t,81] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*g10[t, 6]
    psi[80,t,82] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80,t,83] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80,t,84] <- phi[10,5,t]*psiAa[10,10]*(1-p10[6,t])
    
    
    # Define probabilities of O(t)
    for (j in 1:8){
    for (cc in c(2:4,6)){
    po[((j-1)*6)+cc,t] <- p[j,cc,t] }
    po[1+(6*(j-1)),t] <- 0 #State that can only occur on first capture
    po[5+(6*(j-1)),t] <- 0 #State that can only occur on first capture
    }
    for (jj in c(50:52, 54:56, 58:60, 63:65, 68:70, 72:74, 76:78, 81:83)){
    po[jj,t] <- 1
    }
    po[49,t] <- 0 #State that can only occur on first capture
    po[62,t] <- 0 #State that can only occur on first capture
    po[53,t] <- 0 #Unobservable state
    po[57,t] <- 0 #Unobservable state
    po[61,t] <- 0 #Unobservable state
    po[66,t] <- 0 #Unobservable state
    po[67,t] <- 0 #State that can only occur on first capture
    po[80,t] <- 0 #State that can only occur on first capture
    po[71,t] <- 0 #Unobservable state
    po[75,t] <- 0 #Unobservable state
    po[79,t] <- 0 #Unobservable state
    po[84,t] <- 0 #Unobservable state
    
    # Calculate probability of non-encounter (dq) and reshape the array for the encounter probabilities      
    for (s in 1:ns){
    dp[s,t,s] <- po[s,t]
    dq[s,t,s] <- 1-po[s,t]
    } # s
    
    for (s in 1:(ns-1)){
    for (m in (s+1):ns){
    dp[s,t,m] <- 0
    dq[s,t,m] <- 0
    } # s
    } # m
    for (s in 2:ns){
    for (m in 1:(s-1)){
    dp[s,t,m] <- 0
    dq[s,t,m] <- 0
    } # s
    } # m
    } # t
    
    # Define the multinomial likelihood
    for (t in 1:((n.occasions-1)*ns)){
    marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
    }
    
    # Define the cell probabilities of the multistate m-array   
    # Define matrix U: product of probabilities of state-transition and non-encounter (this is just done because there is no product function for matrix multiplication in JAGS)
    for (t in 1:(n.occasions-2)){
    U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
    for (j in (t+1):(n.occasions-1)){
    U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,t,] %*% dq[,t,]
    }
    }
    U[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ones
    
    # Diagonal
    for (t in 1:(n.occasions-2)){
    pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*% dp[,t,]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% psi[,j,] %*% dp[,j,]
    }
    }
    pr[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- psi[,n.occasions-1,] %*% dp[,n.occasions-1,]
    
    # Below main diagonal
    for (t in 2:(n.occasions-1)){
    for (j in 1:(t-1)){
    pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
    } #j
    } #t
    
    # Last column: probability of non-recapture
    for (t in 1:((n.occasions-1)*ns)){
    pr[t,(n.occasions*ns-(ns-1))] <- 1-sum(pr[t,1:((n.occasions-1)*ns)])
    } #t
    }    ",fill = TRUE)
sink()

# Bundle data
n.unobs <- 8 #Number of unobservable states
n.notobs <- 15 #Number of observable states that dont appear anywhere in capture histories
ns <- length(unique(as.numeric(mat.enc3))) - 1 + n.unobs + n.notobs # calculate the number of states

jags.data <- list(marr = ms.arr, n.occasions = ncol(mat.enc3), rel = rowSums(ms.arr), ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), dista=distmat, ind=in., p.all.years=p.all.years)

inits <- function(){list(a.del=rep(0.5, length.out=3),   
                         mean.phi = 0.1, sigma.phi = runif(1, 0, 5),
                         mean.p = 0.1, sigma.p = runif(1, 0, 5), A.delp=rnorm(3),
                         t.g9=rnorm(4), a.g.3 = 0.5,a.g.4 = 0.5,a.g.6 = 0.5, t.g10 = 0.5,
                         B=-0.1,
                         theta2 = rnorm(1), theta3 = rnorm(1), theta4 = rnorm(1), thetaAa = rnorm(1), thetaAb = rnorm(1), repuls=rnorm(10), attract=rnorm(10)
)}  

# Parameters monitored
parameters <- c("phi.real", "sigma2.phi.real", "p.real", "sigma2.p.real", "p9.real","g9.real", "p10.real", "g10.real","B", "psi1", "psi2","psi3", "psiAa", "psiAb")

# MCMC settings
ni <- 1500 #000 
nt <- 1
nb <- 100 #000 
nc <- 3

# Call JAGS from R, using jagsUI syntax
starttime <- Sys.time()
m9 <- jags(jags.data, inits, parameters, "model9.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=FALSE)
endtime <- Sys.time()
options(max.print=999999) #Default is 99999
print(m9, digits = 3) # B=-0.13  phi for adults still low
print(m9$summary[1:400,], digits = 3)










########################
########################
#Model 10
########################
#Two different mu for p

#Creating an index matrix, to use in script
in. <- matrix(data=rep(1:9,10), ncol=9, nrow=10, byrow=T)
for (k in 1:9){in.[k,k] <- 10}

#Study characteristics
n.occasions <- 10 #Using the data from the last ten years only 
#(thus some birds will appear as they were first marked as 2/3/4-years-olds etc)
n.states <- 85 #Excluding the Elsewhere states
n.obs <- 85 

#Removing the first years of data
enc0817 <- enc[, c(1:5, 19:31)] 

#Combining the numbers for adults as chicks as information about ringing stage is incorporated in the state
enc0817$No. <- enc0817$Ch + enc0817$Ad
#Selecting the relevant columns and changing the order
enc0817. <- enc0817[, c(18, 2:15, 19)] 

#Creating one encounter history for each individual
mat.enc <- matrix(NA, ncol = n.occasions, nrow = sum(enc0817.[,ncol(enc0817.)]))
i <- 1
for (l in 1:enc0817.[i,ncol(enc0817.)]){
  for (k in 1:n.occasions){
    mat.enc[l,k] <- enc0817.[i, (ncol(enc0817.)-n.occasions+k-1)]  } }

for (i in 2:nrow(enc0817.)){
  j <- 1 + sum(enc0817.[1:(i-1),ncol(enc0817.)])
  for (l in j:(sum(enc0817.[1:(i-1),ncol(enc0817.)])+enc0817.[i,ncol(enc0817.)])){
    for (k in 1:n.occasions){
      mat.enc[l,k] <- enc0817.[i, (ncol(enc0817.)-n.occasions+k-1)]    } }}

#Removing the rows with no events (only 0)
mat.enc2 <- mat.enc[which(rowSums(mat.enc)>0),]

#Checking the encounter histories
table(mat.enc2) #Need to check if there are any events that didn't happen, as some probabilities are low

#Recoding the Elsewhere states to 0
elsew <- which(mat.enc2==88|mat.enc2==85|mat.enc2==86)
mat.enc3 <- mat.enc2
mat.enc3[elsew] <- 0

#Checking rows with no events (only 0) 
mat.enc3[which(rowSums(mat.enc3)==0),]
mat.enc2[which(rowSums(mat.enc3)==0),]
#Removing the rows with no events (only 0) - individuals marked prior to first year analyzed, and only encountered Elsewhere, hence with only 0 in encounter history
mat.enc3 <- mat.enc3[which(rowSums(mat.enc3)>0),]
# Compute vector with occasion of first capture
f <- apply(mat.enc3, 1, get.first)

#Checking the encounter histories
table(mat.enc3) 
length(table(mat.enc3))
#Reformating the data using the m-array function. 
ms.arr <- marray(mat.enc3, unobs=23) #Unobs needs to be no. of unobservable states + states not observed (despite observable)

p.all.years <- c(1,3,8) #Sites 1,3,8 should have p estimated for all years of the study

#Index for period for p (depending on ringer, second ringer started in 2013, increased effort from 2014)
in.p <- c(rep.int(1,5), rep.int(2,4))

# Analysis of the model
# Specify model in BUGS language
sink("model10.jags") 
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi[fro, stage, t]: survival probability at age/stage, having been at site fro, in period t
    # psi[fro, till]: movement probability from site fro to site till, only 10 of these as stage at next time-step is known and thereby constraints options to the 10 sites
    # p[...]: encounter probability. 
    # g: proportion of birds encountered (observed, captured or observed & captured) that where captured
    # -------------------------------------------------
    # States (S):
    # 1 alive at a1 s1
    # 2 alive at a2 s1
    # 3 alive at a3 s1
    # 4 alive at a4 s1
    # 5 alive at Ad_recent s1
    # 6 alive at Ad_previous s1
    # 7 alive at a1 s2
    # ...
    # 48 alive at Ad_previous s8
    # 49 alive & captured at a1 s9
    # 50 alive & captured at a2 s9
    # 51 alive & observed at a2 s9
    # 52 alive & observed & captured at a2 s9
    # 53 alive but not seen at a2 s9
    # 54 alive & captured at a3 s9
    # ...
    # 84 alive but not seen at Ad_previous s10
    # 85 dead
    #
    # Observations (O):
    # 1 seen at a1 s1
    # 2 seen at a2 s1
    # 3 seen at a3 s1
    # 4 seen at a4 s1
    # 5 seen at Ad s1
    # 6 seen at a1 s2
    # ...
    # 85/0 not seen
    # Unobservable states: 53,57,61,66, 71,75,79,84
    
    # -------------------------------------------------
    
    # Priors and constraints
    # Survival: random effects of time and island, additive effect of age
    for (fro in 1:10){
    for (stage in 1:6){ 
    for (t in 1:(n.occasions-1)){
    logit(phi[fro, stage, t]) <- mu.phi + eps.phi[fro, t] + del[stage] 
    }}}
    for (fro in 1:10){
    for (t in 1:(n.occasions-1)){
    eps.phi[fro, t] ~ dnorm(0, tau.phi)}}
    mean.phi ~ dunif(0, 1)                    # Prior for mean survival
    mu.phi <- log(mean.phi / (1-mean.phi))    # Logit transformation
    sigma.phi ~ dunif(0, 5)                       # Prior for standard deviation
    tau.phi <- pow(sigma.phi, -2)
    sigma2.phi <- pow(sigma.phi, 2)
    sigma2.phi.real <- sigma2.phi * pow(mean.phi, 2) * pow((1-mean.phi), 2) # Temporal variance on real scale
    
    del[1] <- 0                             # Priors for age/stage-effect
    for (stag in 2:4){ 
    del[stag] <- a.del[1]}
    a.del[1] ~ dnorm(0, 0.01)
    for (stag in 5:6){ 
    del[stag] <- a.del[stag-3]
    a.del[stag-3] ~ dnorm(0, 0.01)}
    for (fron in 1:10){  #Back-transform to the probability scale
    for (stage in 1:6){ 
    for (t in 1:(n.occasions-1)){
    phi.real[fron, stage, t] <- 1/(1+exp(-(mu.phi + eps.phi[fron,t] + del[stage])))    }}}
    
    # Encounters at site 1-8: two intercepts (depending on period - different ringers), random effect of time and island, additive effect of age/stage
    #Sites 1,3,8: 
    for (isl in p.all.years){ 
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[isl, stage, t]) <- mu.p[in.p[t]] + eps.p[isl, t] + delp[stage]     }}}
    for (isl in p.all.years){
    for (t in 1:(n.occasions-1)){
    eps.p[isl, t] ~ dnorm(0, tau.p)}}
    for (isln in p.all.years){  #Back-transform to the probability scale
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    p.real[isln, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[isln,t] + delp[stage])))    }}}
    #Site 2 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[2, stage, t]) <- mu.p[in.p[t]] + eps.p[2, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[2, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[2, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[2,t] + delp[stage])))    }}
    #Site 4 (should have p estimated for 2009-2011,2013-2017):
    for (stage in c(2:4,6)){ 
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    logit(p[4, stage, t]) <- mu.p[in.p[t]] + eps.p[4, t] + delp[stage]   }}
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    eps.p[4, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    p.real[4, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[4,t] + delp[stage])))    }}
    for (stage in c(2:4,6)){ 
    for (t in (n.occasions-6)){
    p[4, stage, t] <- 0   }}
    #Site 5 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[5, stage, t]) <- mu.p[in.p[t]] + eps.p[5, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[5, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[5, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[5,t] + delp[stage])))    }}
    #Site 6 (should have p estimated for 2011-2017):
    for (stage in c(2:4,6)){ 
    for (t in (n.occasions-7):(n.occasions-1)){
    logit(p[6, stage, t]) <- mu.p[in.p[t]] + eps.p[6, t] + delp[stage]   }}
    for (t in (n.occasions-7):(n.occasions-1)){
    eps.p[6, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in (n.occasions-7):(n.occasions-1)){
    p.real[6, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[6,t] + delp[stage])))    }}
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-8)){
    p[6, stage, t] <- 0   }}
    #Site 7 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[7, stage, t]) <- mu.p[in.p[t]] + eps.p[7, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[7, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[7, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[7,t] + delp[stage])))    }}
    
    for (u in 1:2){
    mean.p[u] ~ dunif(0, 1)                       # Prior for mean re-encounter
    mu.p[u] <- log(mean.p[u] / (1-mean.p[u]))}    # Logit transformation
    sigma.p ~ dunif(0, 5)                         # Prior for standard deviation
    tau.p <- pow(sigma.p, -2)
    sigma2.p <- pow(sigma.p, 2)
    sigma2.p.real <- sigma2.p * pow(mean.p, 2) * pow((1-mean.p), 2) # Temporal variance on real scale
    #delp[1] <- 0 #Re-encounters do not happen for age 1 [Shouldn't be necessary to define when using stage c(2:4,6) everywhere for p]
    #delp[5] <- 0 #Re-encounters do not happen for full-grown, just marked [Shouldn't be necessary to define when using stage c(2:4,6) everywhere for p]
    delp[2] <- 0 #Baseline age-class, which mu.p[in.p[t]] applies to
    for (stag in c(3,4)){ # Priors for age/stage-effect
    delp[stag] <- A.delp[stag-2]
    A.delp[stag-2] ~ dnorm(0, 0.01)}
    delp[6] <- A.delp[3]
    A.delp[3] ~ dnorm(0, 0.01)
    
    # Encounters at site 9-10:
    #Site 9 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p9[stage, t]) <- mu.p[in.p[t]] + eps.p9[t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p9[t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p9.real[stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p9[t] + delp[stage])))    }}
    #Site 10 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p10[stage, t]) <- mu.p[in.p[t]] + eps.p10[t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p10[t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p10.real[stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p10[t] + delp[stage])))    }}
    
    for (fd in 1:(n.occasions-5)){
    for (st in c(2:4,6)){
    g9[fd, st] <- 1  }}
    for (fd in (n.occasions-4):(n.occasions-1)){
    logit(g9[fd, 2]) <- t.g9[fd-n.occasions+5]
    logit(g9[fd, 3]) <- t.g9[fd-n.occasions+5] + a.g.3
    logit(g9[fd, 4]) <- t.g9[fd-n.occasions+5] + a.g.4
    logit(g9[fd, 6]) <- t.g9[fd-n.occasions+5] + a.g.6
    t.g9[fd-n.occasions+5] ~ dnorm(0, 0.01)
    g9.real[fd,2] <- 1/(1+exp(-(t.g9[fd-n.occasions+5]))) 
    g9.real[fd,3] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.3)))
    g9.real[fd,4] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.4)))
    g9.real[fd,6] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.6)))}
    for (fd in 1:(n.occasions-4)){ 
    for (st in c(2:4,6)){
    g10[fd, st] <- 1  }}
    logit(g10[n.occasions-3, 2]) <- t.g10
    logit(g10[n.occasions-3, 3]) <- t.g10 + a.g.3
    logit(g10[n.occasions-3, 4]) <- t.g10 + a.g.4
    logit(g10[n.occasions-3, 6]) <- t.g10 + a.g.6
    t.g10 ~ dnorm(0, 0.01)
    g10.real[n.occasions-3,2] <- 1/(1+exp(-(t.g10))) 
    g10.real[n.occasions-3,3] <- 1/(1+exp(-(t.g10 + a.g.3)))
    g10.real[n.occasions-3,4] <- 1/(1+exp(-(t.g10 + a.g.4)))
    g10.real[n.occasions-3,6] <- 1/(1+exp(-(t.g10 + a.g.6)))
    for (fd in (n.occasions-2):(n.occasions-1)){ 
    for (st in c(2:4,6)){
    g10[fd, st] <- 1  }}
    a.g.3 ~ dnorm(0, 0.01)
    a.g.4 ~ dnorm(0, 0.01)
    a.g.6 ~ dnorm(0, 0.01)
    
    # Transitions: multinomial logit
    # Normal priors on logit of all but one transition probs
    for (fro in 1:10){
    for (till in 1:(fro-1)){
    lpsi1[fro, till] <- repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi2[fro, till] <- theta2 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi3[fro, till] <- theta3 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi4[fro, till] <- theta4 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAa[fro, till] <- thetaAa + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAb[fro, till] <- thetaAb + repuls[fro] + attract[till] + B*dista[fro,till]}
    for (till in (fro+1):10){
    lpsi1[fro, till] <- repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi2[fro, till] <- theta2 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi3[fro, till] <- theta3 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi4[fro, till] <- theta4 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAa[fro, till] <- thetaAa + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAb[fro, till] <- thetaAb + repuls[fro] + attract[till] + B*dista[fro,till]} }
    
    for (fro in 1:10){
    for (till in 1:(fro-1)){
    psi1[fro, till] <- exp(lpsi1[fro, till]) / (1 + exp(lpsi1[fro,ind[fro,1]]) + exp(lpsi1[fro,ind[fro,2]]) + exp(lpsi1[fro,ind[fro,3]]) + exp(lpsi1[fro,ind[fro,4]]) + exp(lpsi1[fro,ind[fro,5]]) + exp(lpsi1[fro,ind[fro,6]]) + exp(lpsi1[fro,ind[fro,7]]) + exp(lpsi1[fro,ind[fro,8]]) + exp(lpsi1[fro,ind[fro,9]]))
    psi2[fro, till] <- exp(lpsi2[fro, till]) / (1 + exp(lpsi2[fro,ind[fro,1]]) + exp(lpsi2[fro,ind[fro,2]]) + exp(lpsi2[fro,ind[fro,3]]) + exp(lpsi2[fro,ind[fro,4]]) + exp(lpsi2[fro,ind[fro,5]]) + exp(lpsi2[fro,ind[fro,6]]) + exp(lpsi2[fro,ind[fro,7]]) + exp(lpsi2[fro,ind[fro,8]]) + exp(lpsi2[fro,ind[fro,9]]))
    psi3[fro, till] <- exp(lpsi3[fro, till]) / (1 + exp(lpsi3[fro,ind[fro,1]]) + exp(lpsi3[fro,ind[fro,2]]) + exp(lpsi3[fro,ind[fro,3]]) + exp(lpsi3[fro,ind[fro,4]]) + exp(lpsi3[fro,ind[fro,5]]) + exp(lpsi3[fro,ind[fro,6]]) + exp(lpsi3[fro,ind[fro,7]]) + exp(lpsi3[fro,ind[fro,8]]) + exp(lpsi3[fro,ind[fro,9]]))
    psi4[fro, till] <- exp(lpsi4[fro, till]) / (1 + exp(lpsi4[fro,ind[fro,1]]) + exp(lpsi4[fro,ind[fro,2]]) + exp(lpsi4[fro,ind[fro,3]]) + exp(lpsi4[fro,ind[fro,4]]) + exp(lpsi4[fro,ind[fro,5]]) + exp(lpsi4[fro,ind[fro,6]]) + exp(lpsi4[fro,ind[fro,7]]) + exp(lpsi4[fro,ind[fro,8]]) + exp(lpsi4[fro,ind[fro,9]]))
    psiAa[fro, till] <- exp(lpsiAa[fro, till]) / (1 + exp(lpsiAa[fro,ind[fro,1]]) + exp(lpsiAa[fro,ind[fro,2]]) + exp(lpsiAa[fro,ind[fro,3]]) + exp(lpsiAa[fro,ind[fro,4]]) + exp(lpsiAa[fro,ind[fro,5]]) + exp(lpsiAa[fro,ind[fro,6]]) + exp(lpsiAa[fro,ind[fro,7]]) + exp(lpsiAa[fro,ind[fro,8]]) + exp(lpsiAa[fro,ind[fro,9]]))
    psiAb[fro, till] <- exp(lpsiAb[fro, till]) / (1 + exp(lpsiAb[fro,ind[fro,1]]) + exp(lpsiAb[fro,ind[fro,2]]) + exp(lpsiAb[fro,ind[fro,3]]) + exp(lpsiAb[fro,ind[fro,4]]) + exp(lpsiAb[fro,ind[fro,5]]) + exp(lpsiAb[fro,ind[fro,6]]) + exp(lpsiAb[fro,ind[fro,7]]) + exp(lpsiAb[fro,ind[fro,8]]) + exp(lpsiAb[fro,ind[fro,9]]))  }
    for (till in (fro+1):10){
    psi1[fro, till] <- exp(lpsi1[fro, till]) / (1 + exp(lpsi1[fro,ind[fro,1]]) + exp(lpsi1[fro,ind[fro,2]]) + exp(lpsi1[fro,ind[fro,3]]) + exp(lpsi1[fro,ind[fro,4]]) + exp(lpsi1[fro,ind[fro,5]]) + exp(lpsi1[fro,ind[fro,6]]) + exp(lpsi1[fro,ind[fro,7]]) + exp(lpsi1[fro,ind[fro,8]]) + exp(lpsi1[fro,ind[fro,9]]))
    psi2[fro, till] <- exp(lpsi2[fro, till]) / (1 + exp(lpsi2[fro,ind[fro,1]]) + exp(lpsi2[fro,ind[fro,2]]) + exp(lpsi2[fro,ind[fro,3]]) + exp(lpsi2[fro,ind[fro,4]]) + exp(lpsi2[fro,ind[fro,5]]) + exp(lpsi2[fro,ind[fro,6]]) + exp(lpsi2[fro,ind[fro,7]]) + exp(lpsi2[fro,ind[fro,8]]) + exp(lpsi2[fro,ind[fro,9]]))
    psi3[fro, till] <- exp(lpsi3[fro, till]) / (1 + exp(lpsi3[fro,ind[fro,1]]) + exp(lpsi3[fro,ind[fro,2]]) + exp(lpsi3[fro,ind[fro,3]]) + exp(lpsi3[fro,ind[fro,4]]) + exp(lpsi3[fro,ind[fro,5]]) + exp(lpsi3[fro,ind[fro,6]]) + exp(lpsi3[fro,ind[fro,7]]) + exp(lpsi3[fro,ind[fro,8]]) + exp(lpsi3[fro,ind[fro,9]]))
    psi4[fro, till] <- exp(lpsi4[fro, till]) / (1 + exp(lpsi4[fro,ind[fro,1]]) + exp(lpsi4[fro,ind[fro,2]]) + exp(lpsi4[fro,ind[fro,3]]) + exp(lpsi4[fro,ind[fro,4]]) + exp(lpsi4[fro,ind[fro,5]]) + exp(lpsi4[fro,ind[fro,6]]) + exp(lpsi4[fro,ind[fro,7]]) + exp(lpsi4[fro,ind[fro,8]]) + exp(lpsi4[fro,ind[fro,9]]))
    psiAa[fro, till] <- exp(lpsiAa[fro, till]) / (1 + exp(lpsiAa[fro,ind[fro,1]]) + exp(lpsiAa[fro,ind[fro,2]]) + exp(lpsiAa[fro,ind[fro,3]]) + exp(lpsiAa[fro,ind[fro,4]]) + exp(lpsiAa[fro,ind[fro,5]]) + exp(lpsiAa[fro,ind[fro,6]]) + exp(lpsiAa[fro,ind[fro,7]]) + exp(lpsiAa[fro,ind[fro,8]]) + exp(lpsiAa[fro,ind[fro,9]]))
    psiAb[fro, till] <- exp(lpsiAb[fro, till]) / (1 + exp(lpsiAb[fro,ind[fro,1]]) + exp(lpsiAb[fro,ind[fro,2]]) + exp(lpsiAb[fro,ind[fro,3]]) + exp(lpsiAb[fro,ind[fro,4]]) + exp(lpsiAb[fro,ind[fro,5]]) + exp(lpsiAb[fro,ind[fro,6]]) + exp(lpsiAb[fro,ind[fro,7]]) + exp(lpsiAb[fro,ind[fro,8]]) + exp(lpsiAb[fro,ind[fro,9]]))  }
    #The case below is when fro=till, but rewritten to match JAGS syntax
    psi1[fro,fro] <- 1- sum(psi1[fro,ind[fro,]]) 
    psi2[fro,fro] <- 1- sum(psi2[fro,ind[fro,]]) 
    psi3[fro,fro] <- 1- sum(psi3[fro,ind[fro,]]) 
    psi4[fro,fro] <- 1- sum(psi4[fro,ind[fro,]]) 
    psiAa[fro,fro] <- 1- sum(psiAa[fro,ind[fro,]])   
    psiAb[fro,fro] <- 1- sum(psiAb[fro,ind[fro,]])   
    } 
    
    #Priors for transitions
    theta2 ~ dnorm(0, 0.01)               # Prior for age-effect in psi
    theta3 ~ dnorm(0, 0.01)  
    theta4 ~ dnorm(0, 0.01)  
    thetaAa ~ dnorm(0, 0.01)  
    thetaAb ~ dnorm(0, 0.01) 
    #mean.psi ~ dunif(0, 1)                # Prior for mean movement prob - not used as no transition makes sense to set as baseline
    #mu <- log(mean.psi / (1-mean.psi))    # Logit transformation
    for (w in 1:10){
    repuls[w] ~ dnorm(0, 0.01) 
    attract[w] ~ dnorm(0, 0.01) 
    }
    B ~ dnorm(0, 0.001)I(-10,10)
    
    # Define state-transition and observation matrices 	
    # Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
    for (d in 1:84){
    for (m in 1:8){
    psi[d,t,1+(6*(m-1))] <- 0
    psi[d,t,5+(6*(m-1))] <- 0
    } #m
    psi[d,t,49] <- 0
    psi[d,t,62] <- 0
    psi[d,t,67] <- 0
    psi[d,t,80] <- 0
    } #d
    
    for (sf in 1:8){  
    for (s in 1:8){
    psi[(6*(sf-1))+1,t,(s*6-4)] <- phi[sf,1,t]*psi1[sf,s]
    psi[(6*(sf-1))+1,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+1,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+1,t,(s*6)] <- 0
    psi[(6*(sf-1))+2,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+2,t,(s*6-3)] <- phi[sf,2,t]*psi2[sf,s]
    psi[(6*(sf-1))+2,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+2,t,(s*6)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-2)] <- phi[sf,3,t]*psi3[sf,s]
    psi[(6*(sf-1))+3,t,(s*6)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+4,t,(s*6)] <- phi[sf,4,t]*psi4[sf,s]
    psi[(6*(sf-1))+5,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+5,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+5,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+5,t,(s*6)] <- phi[sf,5,t]*psiAa[sf,s]
    psi[(6*(sf-1))+6,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+6,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+6,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+6,t,(s*6)] <- phi[sf,6,t]*psiAb[sf,s]
    } #s
    
    psi[(6*(sf-1))+1,t,50] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*g9[t, 2]
    psi[(6*(sf-1))+1,t,51] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[(6*(sf-1))+1,t,52] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[(6*(sf-1))+1,t,53] <- phi[sf,1,t]*psi1[sf,9]*(1-p9[2,t])
    for (ex in c(54:61,63:66)){
    psi[(6*(sf-1))+1,t,ex] <- 0  }
    psi[(6*(sf-1))+1,t,68] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*g10[t, 2]
    psi[(6*(sf-1))+1,t,69] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[(6*(sf-1))+1,t,70] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[(6*(sf-1))+1,t,71] <- phi[sf,1,t]*psi1[sf,10]*(1-p10[2,t])
    for (ex in c(72:79,81:84)){
    psi[(6*(sf-1))+1,t,ex] <- 0  }
    
    for (ew in 50:53){
    psi[(6*(sf-1))+2,t,ew] <- 0}
    psi[(6*(sf-1))+2,t,54] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*g9[t, 3]
    psi[(6*(sf-1))+2,t,55] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[(6*(sf-1))+2,t,56] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[(6*(sf-1))+2,t,57] <- phi[sf,2,t]*psi2[sf,9]*(1-p9[3,t])
    for (er in c(58:61,63:66)){
    psi[(6*(sf-1))+2,t,er] <- 0  }
    for (ew in 68:71){
    psi[(6*(sf-1))+2,t,ew] <- 0}
    psi[(6*(sf-1))+2,t,72] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*g10[t, 3]
    psi[(6*(sf-1))+2,t,73] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[(6*(sf-1))+2,t,74] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[(6*(sf-1))+2,t,75] <- phi[sf,2,t]*psi2[sf,10]*(1-p10[3,t])
    for (er in c(76:79,81:84)){
    psi[(6*(sf-1))+2,t,er] <- 0  }
    
    for (ed in 50:57){
    psi[(6*(sf-1))+3,t,ed] <- 0}
    psi[(6*(sf-1))+3,t,58] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*g9[t, 4]
    psi[(6*(sf-1))+3,t,59] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[(6*(sf-1))+3,t,60] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[(6*(sf-1))+3,t,61] <- phi[sf,3,t]*psi3[sf,9]*(1-p9[4,t])
    for (ef in 63:66){
    psi[(6*(sf-1))+3,t,ef] <- 0 }
    for (ed in 68:75){
    psi[(6*(sf-1))+3,t,ed] <- 0}
    psi[(6*(sf-1))+3,t,76] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*g10[t, 4]
    psi[(6*(sf-1))+3,t,77] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[(6*(sf-1))+3,t,78] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[(6*(sf-1))+3,t,79] <- phi[sf,3,t]*psi3[sf,10]*(1-p10[4,t])
    for (ef in 81:84){
    psi[(6*(sf-1))+3,t,ef] <- 0 }
    
    for (es in 50:61){
    psi[(6*(sf-1))+4,t,es] <- 0}
    psi[(6*(sf-1))+4,t,63] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+4,t,64] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+4,t,65] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+4,t,66] <- phi[sf,4,t]*psi4[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+4,t,es] <- 0}
    psi[(6*(sf-1))+4,t,81] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+4,t,82] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+4,t,83] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+4,t,84] <- phi[sf,4,t]*psi4[sf,10]*(1-p10[6,t])
    
    for (es in 50:61){
    psi[(6*(sf-1))+5,t,es] <- 0}
    psi[(6*(sf-1))+5,t,63] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+5,t,64] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+5,t,65] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+5,t,66] <- phi[sf,5,t]*psiAa[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+5,t,es] <- 0}
    psi[(6*(sf-1))+5,t,81] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+5,t,82] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+5,t,83] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+5,t,84] <- phi[sf,5,t]*psiAa[sf,10]*(1-p10[6,t])
    
    for (es in 50:61){
    psi[(6*(sf-1))+6,t,es] <- 0}
    psi[(6*(sf-1))+6,t,63] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+6,t,64] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+6,t,65] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+6,t,66] <- phi[sf,6,t]*psiAb[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+6,t,es] <- 0}
    psi[(6*(sf-1))+6,t,81] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+6,t,82] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+6,t,83] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+6,t,84] <- phi[sf,6,t]*psiAb[sf,10]*(1-p10[6,t])
    } #sf
    
    
    for (s in 1:8){
    psi[49,t,(s*6-4)] <- phi[9,1,t]*psi1[9,s]
    psi[49,t,(s*6-3)] <- 0
    psi[49,t,(s*6-2)] <- 0
    psi[49,t,(s*6)] <- 0
    for (i in 1:4){
    psi[49+i,t,(s*6-4)] <- 0
    psi[49+i,t,(s*6-3)] <- phi[9,2,t]*psi2[9,s]
    psi[49+i,t,(s*6-2)] <- 0
    psi[49+i,t,(s*6)] <- 0
    psi[53+i,t,(s*6-4)] <- 0
    psi[53+i,t,(s*6-3)] <- 0
    psi[53+i,t,(s*6-2)] <- phi[9,3,t]*psi3[9,s]
    psi[53+i,t,(s*6)] <- 0
    psi[57+i,t,(s*6-4)] <- 0
    psi[57+i,t,(s*6-3)] <- 0
    psi[57+i,t,(s*6-2)] <- 0
    psi[57+i,t,(s*6)] <- phi[9,4,t]*psi4[9,s] 
    psi[62+i,t,(s*6-4)] <- 0
    psi[62+i,t,(s*6-3)] <- 0
    psi[62+i,t,(s*6-2)] <- 0
    psi[62+i,t,(s*6)] <- phi[9,6,t]*psiAb[9,s]
    } #i
    psi[62,t,(s*6-4)] <- 0
    psi[62,t,(s*6-3)] <- 0
    psi[62,t,(s*6-2)] <- 0
    psi[62,t,(s*6)] <- phi[9,5,t]*psiAa[9,s]
    } #s
    
    for (s in 1:8){
    psi[67,t,(s*6-4)] <- phi[10,1,t]*psi1[10,s]
    psi[67,t,(s*6-3)] <- 0
    psi[67,t,(s*6-2)] <- 0
    psi[67,t,(s*6)] <- 0
    for (i in 1:4){
    psi[67+i,t,(s*6-4)] <- 0
    psi[67+i,t,(s*6-3)] <- phi[10,2,t]*psi2[10,s]
    psi[67+i,t,(s*6-2)] <- 0
    psi[67+i,t,(s*6)] <- 0
    psi[71+i,t,(s*6-4)] <- 0
    psi[71+i,t,(s*6-3)] <- 0
    psi[71+i,t,(s*6-2)] <- phi[10,3,t]*psi3[10,s]
    psi[71+i,t,(s*6)] <- 0
    psi[75+i,t,(s*6-4)] <- 0
    psi[75+i,t,(s*6-3)] <- 0
    psi[75+i,t,(s*6-2)] <- 0
    psi[75+i,t,(s*6)] <- phi[10,4,t]*psi4[10,s] 
    psi[80+i,t,(s*6-4)] <- 0
    psi[80+i,t,(s*6-3)] <- 0
    psi[80+i,t,(s*6-2)] <- 0
    psi[80+i,t,(s*6)] <- phi[10,6,t]*psiAb[10,s]
    } #i
    psi[80,t,(s*6-4)] <- 0
    psi[80,t,(s*6-3)] <- 0
    psi[80,t,(s*6-2)] <- 0
    psi[80,t,(s*6)] <- phi[10,5,t]*psiAa[10,s]
    } #s
    
    psi[49,t,50] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*g9[t, 2]
    psi[49,t,51] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[49,t,52] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[49,t,53] <- phi[9,1,t]*psi1[9,9]*(1-p9[2,t])
    for (ex in c(54:61,63:66)){
    psi[49,t,ex] <- 0  }
    
    for (i in 1:4){
    for (ew in 50:53){
    psi[49+i,t,ew] <- 0  }
    psi[49+i,t,54] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*g9[t, 3]
    psi[49+i,t,55] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[49+i,t,56] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[49+i,t,57] <- phi[9,2,t]*psi2[9,9]*(1-p9[3,t])
    for (er in c(58:61,63:66)){
    psi[49+i,t,er] <- 0  }
    for (ed in 50:57){
    psi[53+i,t,ed] <- 0}
    psi[53+i,t,58] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*g9[t, 4]
    psi[53+i,t,59] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[53+i,t,60] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[53+i,t,61] <- phi[9,3,t]*psi3[9,9]*(1-p9[4,t])
    for (ef in 63:66){
    psi[53+i,t,ef] <- 0 }
    for (es in 50:61){
    psi[57+i,t,es] <- 0}
    psi[57+i,t,63] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*g9[t, 6]
    psi[57+i,t,64] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[57+i,t,65] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[57+i,t,66] <- phi[9,4,t]*psi4[9,9]*(1-p9[6,t])
    for (es in 50:61){
    psi[62+i,t,es] <- 0}
    psi[62+i,t,63] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*g9[t, 6]
    psi[62+i,t,64] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62+i,t,65] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62+i,t,66] <- phi[9,6,t]*psiAb[9,9]*(1-p9[6,t])
    } #i
    
    for (es in 50:61){
    psi[62,t,es] <- 0}
    psi[62,t,63] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*g9[t, 6]
    psi[62,t,64] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62,t,65] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62,t,66] <- phi[9,5,t]*psiAa[9,9]*(1-p9[6,t])
    
    
    psi[49,t,50+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*g10[t, 2]
    psi[49,t,51+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[49,t,52+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[49,t,53+18] <- phi[9,1,t]*psi1[9,10]*(1-p10[2,t])
    for (ex in c(54:61,63:66)){
    psi[49,t,ex+18] <- 0  }
    
    for (i in 1:4){
    for (ew in 50:53){
    psi[49+i,t,ew+18] <- 0  }
    psi[49+i,t,54+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*g10[t, 3]
    psi[49+i,t,55+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[49+i,t,56+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[49+i,t,57+18] <- phi[9,2,t]*psi2[9,10]*(1-p10[3,t])
    for (er in c(58:61,63:66)){
    psi[49+i,t,er+18] <- 0  }
    for (ed in 50:57){
    psi[53+i,t,ed+18] <- 0}
    psi[53+i,t,58+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*g10[t, 4]
    psi[53+i,t,59+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[53+i,t,60+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[53+i,t,61+18] <- phi[9,3,t]*psi3[9,10]*(1-p10[4,t])
    for (ef in 63:66){
    psi[53+i,t,ef+18] <- 0 }
    for (es in 50:61){
    psi[57+i,t,es+18] <- 0}
    psi[57+i,t,63+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*g10[t, 6]
    psi[57+i,t,64+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[57+i,t,65+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[57+i,t,66+18] <- phi[9,4,t]*psi4[9,10]*(1-p10[6,t])
    for (es in 50:61){
    psi[62+i,t,es+18] <- 0}
    psi[62+i,t,63+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*g10[t, 6]
    psi[62+i,t,64+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62+i,t,65+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62+i,t,66+18] <- phi[9,6,t]*psiAb[9,10]*(1-p10[6,t])
    } #i
    
    for (es in 50:61){
    psi[62,t,es+18] <- 0}
    psi[62,t,63+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*g10[t, 6]
    psi[62,t,64+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62,t,65+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62,t,66+18] <- phi[9,5,t]*psiAa[9,10]*(1-p10[6,t])
    
    psi[67,t,50] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*g9[t, 2]
    psi[67,t,51] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[67,t,52] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[67,t,53] <- phi[10,1,t]*psi1[10,9]*(1-p9[2,t])
    for (ex in c(72:79,81:84)){
    psi[67,t,ex-18] <- 0  }
    
    for (i in 1:4){
    for (ew in 68:71){
    psi[67+i,t,ew-18] <- 0  }
    psi[67+i,t,72-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*g9[t, 3]
    psi[67+i,t,73-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[67+i,t,74-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[67+i,t,75-18] <- phi[10,2,t]*psi2[10,9]*(1-p9[3,t])
    for (er in c(76:79,81:84)){
    psi[67+i,t,er-18] <- 0  }
    for (ed in 68:75){
    psi[71+i,t,ed-18] <- 0}
    psi[71+i,t,76-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*g9[t, 4]
    psi[71+i,t,77-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[71+i,t,78-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[71+i,t,79-18] <- phi[10,3,t]*psi3[10,9]*(1-p9[4,t])
    for (ef in 81:84){
    psi[71+i,t,ef-18] <- 0 }
    for (es in 68:79){
    psi[75+i,t,es-18] <- 0}
    psi[75+i,t,81-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*g9[t, 6]
    psi[75+i,t,82-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[75+i,t,83-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[75+i,t,84-18] <- phi[10,4,t]*psi4[10,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[80+i,t,es-18] <- 0}
    psi[80+i,t,81-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*g9[t, 6]
    psi[80+i,t,82-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80+i,t,83-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80+i,t,84-18] <- phi[10,6,t]*psiAb[10,9]*(1-p9[6,t])
    } #i
    
    for (es in 68:79){
    psi[80,t,es-18] <- 0}
    psi[80,t,81-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*g9[t, 6]
    psi[80,t,82-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80,t,83-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80,t,84-18] <- phi[10,5,t]*psiAa[10,9]*(1-p9[6,t])
    
    
    psi[67,t,68] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*g10[t, 2]
    psi[67,t,69] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[67,t,70] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[67,t,71] <- phi[10,1,t]*psi1[10,10]*(1-p10[2,t])
    for (ex in c(72:79,81:84)){
    psi[67,t,ex] <- 0  }
    
    for (i in 1:4){
    for (ew in 68:71){
    psi[67+i,t,ew] <- 0  }
    psi[67+i,t,72] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*g10[t, 3]
    psi[67+i,t,73] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[67+i,t,74] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[67+i,t,75] <- phi[10,2,t]*psi2[10,10]*(1-p10[3,t])
    for (er in c(76:79,81:84)){
    psi[67+i,t,er] <- 0  }
    for (ed in 68:75){
    psi[71+i,t,ed] <- 0}
    psi[71+i,t,76] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*g10[t, 4]
    psi[71+i,t,77] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[71+i,t,78] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[71+i,t,79] <- phi[10,3,t]*psi3[10,10]*(1-p10[4,t])
    for (ef in 81:84){
    psi[71+i,t,ef] <- 0 }
    for (es in 68:79){
    psi[75+i,t,es] <- 0}
    psi[75+i,t,81] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*g10[t, 6]
    psi[75+i,t,82] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[75+i,t,83] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[75+i,t,84] <- phi[10,4,t]*psi4[10,10]*(1-p10[6,t])
    for (es in 68:79){
    psi[80+i,t,es] <- 0}
    psi[80+i,t,81] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*g10[t, 6]
    psi[80+i,t,82] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80+i,t,83] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80+i,t,84] <- phi[10,6,t]*psiAb[10,10]*(1-p10[6,t])
    } #i
    for (es in 68:79){
    psi[80,t,es] <- 0}
    psi[80,t,81] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*g10[t, 6]
    psi[80,t,82] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80,t,83] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80,t,84] <- phi[10,5,t]*psiAa[10,10]*(1-p10[6,t])
    
    
    # Define probabilities of O(t)
    for (j in 1:8){
    for (cc in c(2:4,6)){
    po[((j-1)*6)+cc,t] <- p[j,cc,t] }
    po[1+(6*(j-1)),t] <- 0 #State that can only occur on first capture
    po[5+(6*(j-1)),t] <- 0 #State that can only occur on first capture
    }
    for (jj in c(50:52, 54:56, 58:60, 63:65, 68:70, 72:74, 76:78, 81:83)){
    po[jj,t] <- 1
    }
    po[49,t] <- 0 #State that can only occur on first capture
    po[62,t] <- 0 #State that can only occur on first capture
    po[53,t] <- 0 #Unobservable state
    po[57,t] <- 0 #Unobservable state
    po[61,t] <- 0 #Unobservable state
    po[66,t] <- 0 #Unobservable state
    po[67,t] <- 0 #State that can only occur on first capture
    po[80,t] <- 0 #State that can only occur on first capture
    po[71,t] <- 0 #Unobservable state
    po[75,t] <- 0 #Unobservable state
    po[79,t] <- 0 #Unobservable state
    po[84,t] <- 0 #Unobservable state
    
    # Calculate probability of non-encounter (dq) and reshape the array for the encounter probabilities      
    for (s in 1:ns){
    dp[s,t,s] <- po[s,t]
    dq[s,t,s] <- 1-po[s,t]
    } # s
    
    for (s in 1:(ns-1)){
    for (m in (s+1):ns){
    dp[s,t,m] <- 0
    dq[s,t,m] <- 0
    } # s
    } # m
    for (s in 2:ns){
    for (m in 1:(s-1)){
    dp[s,t,m] <- 0
    dq[s,t,m] <- 0
    } # s
    } # m
    } # t
    
    # Define the multinomial likelihood
    for (t in 1:((n.occasions-1)*ns)){
    marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
    }
    
    # Define the cell probabilities of the multistate m-array   
    # Define matrix U: product of probabilities of state-transition and non-encounter (this is just done because there is no product function for matrix multiplication in JAGS)
    for (t in 1:(n.occasions-2)){
    U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
    for (j in (t+1):(n.occasions-1)){
    U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,t,] %*% dq[,t,]
    }
    }
    U[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ones
    
    # Diagonal
    for (t in 1:(n.occasions-2)){
    pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*% dp[,t,]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% psi[,j,] %*% dp[,j,]
    }
    }
    pr[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- psi[,n.occasions-1,] %*% dp[,n.occasions-1,]
    
    # Below main diagonal
    for (t in 2:(n.occasions-1)){
    for (j in 1:(t-1)){
    pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
    } #j
    } #t
    
    # Last column: probability of non-recapture
    for (t in 1:((n.occasions-1)*ns)){
    pr[t,(n.occasions*ns-(ns-1))] <- 1-sum(pr[t,1:((n.occasions-1)*ns)])
    } #t
    }    ",fill = TRUE)
sink()

# Bundle data
n.unobs <- 8 #Number of unobservable states
n.notobs <- 15 #Number of observable states that dont appear anywhere in capture histories
ns <- length(unique(as.numeric(mat.enc3))) - 1 + n.unobs + n.notobs # calculate the number of states

jags.data <- list(marr = ms.arr, n.occasions = ncol(mat.enc3), rel = rowSums(ms.arr), ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), 
                  dista=distmat, ind=in., in.p=in.p, p.all.years=p.all.years)

inits <- function(){list(a.del=rep(0.5, length.out=3),   
                         mean.phi = 0.1, sigma.phi = runif(1, 0, 5),
                         mean.p = runif(2,min=0,max=1), sigma.p = runif(1, 0, 5), A.delp=rnorm(3),
                         t.g9=rnorm(4), a.g.3 = 0.5,a.g.4 = 0.5,a.g.6 = 0.5, t.g10 = 0.5,
                         B=rnorm(n=1, mean=-0.3, sd=0.1),
                         theta2 = rnorm(1), theta3 = rnorm(1), theta4 = rnorm(1), thetaAa = rnorm(1), thetaAb = rnorm(1), repuls=rnorm(10), attract=rnorm(10)
)}  

# Parameters monitored
parameters <- c("phi.real", "sigma2.phi.real", "p.real", "sigma2.p.real", "p9.real","g9.real", "p10.real", "g10.real","B", "psi1", "psi2","psi3", "psiAa", "psiAb")

# MCMC settings
ni <- 1500 #0 
nt <- 1
nb <- 200 #0
nc <- 3

# Call JAGS from R, using jagsUI syntax
starttime <- Sys.time()
m10 <- jags(jags.data, inits, parameters, "model10.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
endtime <- Sys.time()
options(max.print=999999) #Default is 99999
print(m10, digits = 3) # 
print(m10$summary[1:400,], digits = 3)









########################
########################
#Model 11
########################
#Different p at site 9 and 10 when resight occur, in addition to recaptures


# Analysis of the model
# Specify model in BUGS language
sink("model11.jags") 
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi[fro, stage, t]: survival probability at age/stage, having been at site fro, in period t
    # psi[fro, till]: movement probability from site fro to site till, only 10 of these as stage at next time-step is known and thereby constraints options to the 10 sites
    # p[...]: encounter probability. 
    # g: proportion of birds encountered (observed, captured or observed & captured) that where captured
    # -------------------------------------------------
    # States (S):
    # 1 alive at a1 s1
    # 2 alive at a2 s1
    # 3 alive at a3 s1
    # 4 alive at a4 s1
    # 5 alive at Ad_recent s1
    # 6 alive at Ad_previous s1
    # 7 alive at a1 s2
    # ...
    # 48 alive at Ad_previous s8
    # 49 alive & captured at a1 s9
    # 50 alive & captured at a2 s9
    # 51 alive & observed at a2 s9
    # 52 alive & observed & captured at a2 s9
    # 53 alive but not seen at a2 s9
    # 54 alive & captured at a3 s9
    # ...
    # 84 alive but not seen at Ad_previous s10
    # 85 dead
    #
    # Observations (O):
    # 1 seen at a1 s1
    # 2 seen at a2 s1
    # 3 seen at a3 s1
    # 4 seen at a4 s1
    # 5 seen at Ad s1
    # 6 seen at a1 s2
    # ...
    # 85/0 not seen
    # Unobservable states: 53,57,61,66, 71,75,79,84
    
    # -------------------------------------------------
    
    # Priors and constraints
    # Survival: random effects of time and island, additive effect of age
    for (fro in 1:10){
    for (stage in 1:6){ 
    for (t in 1:(n.occasions-1)){
    logit(phi[fro, stage, t]) <- mu.phi + eps.phi[fro, t] + del[stage] 
    }}}
    for (fro in 1:10){
    for (t in 1:(n.occasions-1)){
    eps.phi[fro, t] ~ dnorm(0, tau.phi)}}
    mean.phi ~ dunif(0, 1)                    # Prior for mean survival
    mu.phi <- log(mean.phi / (1-mean.phi))    # Logit transformation
    sigma.phi ~ dunif(0, 5)                       # Prior for standard deviation
    tau.phi <- pow(sigma.phi, -2)
    sigma2.phi <- pow(sigma.phi, 2)
    sigma2.phi.real <- sigma2.phi * pow(mean.phi, 2) * pow((1-mean.phi), 2) # Temporal variance on real scale
    
    del[1] <- 0                             # Priors for age/stage-effect
    for (stag in 2:4){ 
    del[stag] <- a.del[1]}
    a.del[1] ~ dnorm(0, 0.01)
    for (stag in 5:6){ 
    del[stag] <- a.del[stag-3]
    a.del[stag-3] ~ dnorm(0, 0.01)}
    for (fron in 1:10){  #Back-transform to the probability scale
    for (stage in 1:6){ 
    for (t in 1:(n.occasions-1)){
    phi.real[fron, stage, t] <- 1/(1+exp(-(mu.phi + eps.phi[fron,t] + del[stage])))    }}}
    
    # Encounters at site 1-8: two intercepts (depending on period - different ringers), random effect of time and island, additive effect of age/stage
    #Sites 1,3,8: 
    for (isl in p.all.years){ 
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[isl, stage, t]) <- mu.p[in.p[t]] + eps.p[isl, t] + delp[stage]     }}}
    for (isl in p.all.years){
    for (t in 1:(n.occasions-1)){
    eps.p[isl, t] ~ dnorm(0, tau.p)}}
    for (isln in p.all.years){  #Back-transform to the probability scale
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    p.real[isln, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[isln,t] + delp[stage])))    }}}
    #Site 2 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[2, stage, t]) <- mu.p[in.p[t]] + eps.p[2, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[2, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[2, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[2,t] + delp[stage])))    }}
    #Site 4 (should have p estimated for 2009-2011,2013-2017):
    for (stage in c(2:4,6)){ 
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    logit(p[4, stage, t]) <- mu.p[in.p[t]] + eps.p[4, t] + delp[stage]   }}
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    eps.p[4, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in c((n.occasions-9):(n.occasions-7), (n.occasions-5):(n.occasions-1))){
    p.real[4, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[4,t] + delp[stage])))    }}
    for (stage in c(2:4,6)){ 
    for (t in (n.occasions-6)){
    p[4, stage, t] <- 0   }}
    #Site 5 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[5, stage, t]) <- mu.p[in.p[t]] + eps.p[5, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[5, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[5, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[5,t] + delp[stage])))    }}
    #Site 6 (should have p estimated for 2011-2017):
    for (stage in c(2:4,6)){ 
    for (t in (n.occasions-7):(n.occasions-1)){
    logit(p[6, stage, t]) <- mu.p[in.p[t]] + eps.p[6, t] + delp[stage]   }}
    for (t in (n.occasions-7):(n.occasions-1)){
    eps.p[6, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in (n.occasions-7):(n.occasions-1)){
    p.real[6, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[6,t] + delp[stage])))    }}
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-8)){
    p[6, stage, t] <- 0   }}
    #Site 7 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p[7, stage, t]) <- mu.p[in.p[t]] + eps.p[7, t] + delp[stage]   }}
    for (t in 1:(n.occasions-1)){
    eps.p[7, t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p.real[7, stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p[7,t] + delp[stage])))    }}
    
    for (u in 1:2){
    mean.p[u] ~ dunif(0, 1)                       # Prior for mean re-encounter
    mu.p[u] <- log(mean.p[u] / (1-mean.p[u]))}    # Logit transformation
    sigma.p ~ dunif(0, 5)                         # Prior for standard deviation
    tau.p <- pow(sigma.p, -2)
    sigma2.p <- pow(sigma.p, 2)
    sigma2.p.real <- sigma2.p * pow(mean.p, 2) * pow((1-mean.p), 2) # Temporal variance on real scale
    #delp[1] <- 0 #Re-encounters do not happen for age 1 [Shouldn't be necessary to define when using stage c(2:4,6) everywhere for p]
    #delp[5] <- 0 #Re-encounters do not happen for full-grown, just marked [Shouldn't be necessary to define when using stage c(2:4,6) everywhere for p]
    delp[2] <- 0 #Baseline age-class, which mu.p[in.p[t]] applies to
    for (stag in c(3,4)){ # Priors for age/stage-effect
    delp[stag] <- A.delp[stag-2]
    A.delp[stag-2] ~ dnorm(0, 0.01)}
    delp[6] <- A.delp[3]
    A.delp[3] ~ dnorm(0, 0.01)
    
    # Encounters at site 9-10:
    #Site 9 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p9[stage, t]) <- mu.p[in.p[t]] + eps.p9[t] + delp[stage]   }}
    for (t in 1:(n.occasions-5)){
    eps.p9[t] ~ dnorm(0, tau.p)}
    for (t in (n.occasions-4):(n.occasions-1)){
    eps.p9[t] <- p.both.enc9[t-n.occasions+5]
    p.both.enc9[t-n.occasions+5] ~ dnorm(0, 0.01)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p9.real[stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p9[t] + delp[stage])))    }}
    #Site 10 (should have p estimated for all years in period 2008-2017):
    for (stage in c(2:4,6)){ 
    for (t in 1:(n.occasions-1)){
    logit(p10[stage, t]) <- mu.p[in.p[t]] + eps.p10[t] + delp[stage]   }}
    for (t in 1:(n.occasions-4)){
    eps.p10[t] ~ dnorm(0, tau.p)}
    eps.p10[(n.occasions-3)] <- p.both.enc10
    p.both.enc10 ~ dnorm(0, 0.01)
    for (t in (n.occasions-2):(n.occasions-1)){
    eps.p10[t] ~ dnorm(0, tau.p)}
    for (stage in c(2:4,6)){ #Back-transform to the probability scale
    for (t in 1:(n.occasions-1)){
    p10.real[stage, t] <- 1/(1+exp(-(mu.p[in.p[t]] + eps.p10[t] + delp[stage])))    }}
    
    for (fd in 1:(n.occasions-5)){
    for (st in c(2:4,6)){
    g9[fd, st] <- 1  }}
    for (fd in (n.occasions-4):(n.occasions-1)){
    logit(g9[fd, 2]) <- t.g9[fd-n.occasions+5]
    logit(g9[fd, 3]) <- t.g9[fd-n.occasions+5] + a.g.3
    logit(g9[fd, 4]) <- t.g9[fd-n.occasions+5] + a.g.4
    logit(g9[fd, 6]) <- t.g9[fd-n.occasions+5] + a.g.6
    t.g9[fd-n.occasions+5] ~ dnorm(0, 0.01)
    g9.real[fd,2] <- 1/(1+exp(-(t.g9[fd-n.occasions+5]))) 
    g9.real[fd,3] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.3)))
    g9.real[fd,4] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.4)))
    g9.real[fd,6] <- 1/(1+exp(-(t.g9[fd-n.occasions+5] + a.g.6)))}
    for (fd in 1:(n.occasions-4)){ 
    for (st in c(2:4,6)){
    g10[fd, st] <- 1  }}
    logit(g10[n.occasions-3, 2]) <- t.g10
    logit(g10[n.occasions-3, 3]) <- t.g10 + a.g.3
    logit(g10[n.occasions-3, 4]) <- t.g10 + a.g.4
    logit(g10[n.occasions-3, 6]) <- t.g10 + a.g.6
    t.g10 ~ dnorm(0, 0.01)
    g10.real[n.occasions-3,2] <- 1/(1+exp(-(t.g10))) 
    g10.real[n.occasions-3,3] <- 1/(1+exp(-(t.g10 + a.g.3)))
    g10.real[n.occasions-3,4] <- 1/(1+exp(-(t.g10 + a.g.4)))
    g10.real[n.occasions-3,6] <- 1/(1+exp(-(t.g10 + a.g.6)))
    for (fd in (n.occasions-2):(n.occasions-1)){ 
    for (st in c(2:4,6)){
    g10[fd, st] <- 1  }}
    a.g.3 ~ dnorm(0, 0.01)
    a.g.4 ~ dnorm(0, 0.01)
    a.g.6 ~ dnorm(0, 0.01)
    
    # Transitions: multinomial logit
    # Normal priors on logit of all but one transition probs
    for (fro in 1:10){
    for (till in 1:(fro-1)){
    lpsi1[fro, till] <- repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi2[fro, till] <- theta2 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi3[fro, till] <- theta3 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi4[fro, till] <- theta4 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAa[fro, till] <- thetaAa + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAb[fro, till] <- thetaAb + repuls[fro] + attract[till] + B*dista[fro,till]}
    for (till in (fro+1):10){
    lpsi1[fro, till] <- repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi2[fro, till] <- theta2 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi3[fro, till] <- theta3 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsi4[fro, till] <- theta4 + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAa[fro, till] <- thetaAa + repuls[fro] + attract[till] + B*dista[fro,till]
    lpsiAb[fro, till] <- thetaAb + repuls[fro] + attract[till] + B*dista[fro,till]} }
    
    for (fro in 1:10){
    for (till in 1:(fro-1)){
    psi1[fro, till] <- exp(lpsi1[fro, till]) / (1 + exp(lpsi1[fro,ind[fro,1]]) + exp(lpsi1[fro,ind[fro,2]]) + exp(lpsi1[fro,ind[fro,3]]) + exp(lpsi1[fro,ind[fro,4]]) + exp(lpsi1[fro,ind[fro,5]]) + exp(lpsi1[fro,ind[fro,6]]) + exp(lpsi1[fro,ind[fro,7]]) + exp(lpsi1[fro,ind[fro,8]]) + exp(lpsi1[fro,ind[fro,9]]))
    psi2[fro, till] <- exp(lpsi2[fro, till]) / (1 + exp(lpsi2[fro,ind[fro,1]]) + exp(lpsi2[fro,ind[fro,2]]) + exp(lpsi2[fro,ind[fro,3]]) + exp(lpsi2[fro,ind[fro,4]]) + exp(lpsi2[fro,ind[fro,5]]) + exp(lpsi2[fro,ind[fro,6]]) + exp(lpsi2[fro,ind[fro,7]]) + exp(lpsi2[fro,ind[fro,8]]) + exp(lpsi2[fro,ind[fro,9]]))
    psi3[fro, till] <- exp(lpsi3[fro, till]) / (1 + exp(lpsi3[fro,ind[fro,1]]) + exp(lpsi3[fro,ind[fro,2]]) + exp(lpsi3[fro,ind[fro,3]]) + exp(lpsi3[fro,ind[fro,4]]) + exp(lpsi3[fro,ind[fro,5]]) + exp(lpsi3[fro,ind[fro,6]]) + exp(lpsi3[fro,ind[fro,7]]) + exp(lpsi3[fro,ind[fro,8]]) + exp(lpsi3[fro,ind[fro,9]]))
    psi4[fro, till] <- exp(lpsi4[fro, till]) / (1 + exp(lpsi4[fro,ind[fro,1]]) + exp(lpsi4[fro,ind[fro,2]]) + exp(lpsi4[fro,ind[fro,3]]) + exp(lpsi4[fro,ind[fro,4]]) + exp(lpsi4[fro,ind[fro,5]]) + exp(lpsi4[fro,ind[fro,6]]) + exp(lpsi4[fro,ind[fro,7]]) + exp(lpsi4[fro,ind[fro,8]]) + exp(lpsi4[fro,ind[fro,9]]))
    psiAa[fro, till] <- exp(lpsiAa[fro, till]) / (1 + exp(lpsiAa[fro,ind[fro,1]]) + exp(lpsiAa[fro,ind[fro,2]]) + exp(lpsiAa[fro,ind[fro,3]]) + exp(lpsiAa[fro,ind[fro,4]]) + exp(lpsiAa[fro,ind[fro,5]]) + exp(lpsiAa[fro,ind[fro,6]]) + exp(lpsiAa[fro,ind[fro,7]]) + exp(lpsiAa[fro,ind[fro,8]]) + exp(lpsiAa[fro,ind[fro,9]]))
    psiAb[fro, till] <- exp(lpsiAb[fro, till]) / (1 + exp(lpsiAb[fro,ind[fro,1]]) + exp(lpsiAb[fro,ind[fro,2]]) + exp(lpsiAb[fro,ind[fro,3]]) + exp(lpsiAb[fro,ind[fro,4]]) + exp(lpsiAb[fro,ind[fro,5]]) + exp(lpsiAb[fro,ind[fro,6]]) + exp(lpsiAb[fro,ind[fro,7]]) + exp(lpsiAb[fro,ind[fro,8]]) + exp(lpsiAb[fro,ind[fro,9]]))  }
    for (till in (fro+1):10){
    psi1[fro, till] <- exp(lpsi1[fro, till]) / (1 + exp(lpsi1[fro,ind[fro,1]]) + exp(lpsi1[fro,ind[fro,2]]) + exp(lpsi1[fro,ind[fro,3]]) + exp(lpsi1[fro,ind[fro,4]]) + exp(lpsi1[fro,ind[fro,5]]) + exp(lpsi1[fro,ind[fro,6]]) + exp(lpsi1[fro,ind[fro,7]]) + exp(lpsi1[fro,ind[fro,8]]) + exp(lpsi1[fro,ind[fro,9]]))
    psi2[fro, till] <- exp(lpsi2[fro, till]) / (1 + exp(lpsi2[fro,ind[fro,1]]) + exp(lpsi2[fro,ind[fro,2]]) + exp(lpsi2[fro,ind[fro,3]]) + exp(lpsi2[fro,ind[fro,4]]) + exp(lpsi2[fro,ind[fro,5]]) + exp(lpsi2[fro,ind[fro,6]]) + exp(lpsi2[fro,ind[fro,7]]) + exp(lpsi2[fro,ind[fro,8]]) + exp(lpsi2[fro,ind[fro,9]]))
    psi3[fro, till] <- exp(lpsi3[fro, till]) / (1 + exp(lpsi3[fro,ind[fro,1]]) + exp(lpsi3[fro,ind[fro,2]]) + exp(lpsi3[fro,ind[fro,3]]) + exp(lpsi3[fro,ind[fro,4]]) + exp(lpsi3[fro,ind[fro,5]]) + exp(lpsi3[fro,ind[fro,6]]) + exp(lpsi3[fro,ind[fro,7]]) + exp(lpsi3[fro,ind[fro,8]]) + exp(lpsi3[fro,ind[fro,9]]))
    psi4[fro, till] <- exp(lpsi4[fro, till]) / (1 + exp(lpsi4[fro,ind[fro,1]]) + exp(lpsi4[fro,ind[fro,2]]) + exp(lpsi4[fro,ind[fro,3]]) + exp(lpsi4[fro,ind[fro,4]]) + exp(lpsi4[fro,ind[fro,5]]) + exp(lpsi4[fro,ind[fro,6]]) + exp(lpsi4[fro,ind[fro,7]]) + exp(lpsi4[fro,ind[fro,8]]) + exp(lpsi4[fro,ind[fro,9]]))
    psiAa[fro, till] <- exp(lpsiAa[fro, till]) / (1 + exp(lpsiAa[fro,ind[fro,1]]) + exp(lpsiAa[fro,ind[fro,2]]) + exp(lpsiAa[fro,ind[fro,3]]) + exp(lpsiAa[fro,ind[fro,4]]) + exp(lpsiAa[fro,ind[fro,5]]) + exp(lpsiAa[fro,ind[fro,6]]) + exp(lpsiAa[fro,ind[fro,7]]) + exp(lpsiAa[fro,ind[fro,8]]) + exp(lpsiAa[fro,ind[fro,9]]))
    psiAb[fro, till] <- exp(lpsiAb[fro, till]) / (1 + exp(lpsiAb[fro,ind[fro,1]]) + exp(lpsiAb[fro,ind[fro,2]]) + exp(lpsiAb[fro,ind[fro,3]]) + exp(lpsiAb[fro,ind[fro,4]]) + exp(lpsiAb[fro,ind[fro,5]]) + exp(lpsiAb[fro,ind[fro,6]]) + exp(lpsiAb[fro,ind[fro,7]]) + exp(lpsiAb[fro,ind[fro,8]]) + exp(lpsiAb[fro,ind[fro,9]]))  }
    #The case below is when fro=till, but rewritten to match JAGS syntax
    psi1[fro,fro] <- 1- sum(psi1[fro,ind[fro,]]) 
    psi2[fro,fro] <- 1- sum(psi2[fro,ind[fro,]]) 
    psi3[fro,fro] <- 1- sum(psi3[fro,ind[fro,]]) 
    psi4[fro,fro] <- 1- sum(psi4[fro,ind[fro,]]) 
    psiAa[fro,fro] <- 1- sum(psiAa[fro,ind[fro,]])   
    psiAb[fro,fro] <- 1- sum(psiAb[fro,ind[fro,]])   
    } 
    
    #Priors for transitions
    theta2 ~ dnorm(0, 0.01)               # Prior for age-effect in psi
    theta3 ~ dnorm(0, 0.01)  
    theta4 ~ dnorm(0, 0.01)  
    thetaAa ~ dnorm(0, 0.01)  
    thetaAb ~ dnorm(0, 0.01) 
    #mean.psi ~ dunif(0, 1)                # Prior for mean movement prob - not used as no transition makes sense to set as baseline
    #mu <- log(mean.psi / (1-mean.psi))    # Logit transformation
    for (w in 1:10){
    repuls[w] ~ dnorm(0, 0.01) 
    attract[w] ~ dnorm(0, 0.01) 
    }
    B ~ dnorm(0, 0.001)I(-10,10)
    
    # Define state-transition and observation matrices 	
    # Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
    for (d in 1:84){
    for (m in 1:8){
    psi[d,t,1+(6*(m-1))] <- 0
    psi[d,t,5+(6*(m-1))] <- 0
    } #m
    psi[d,t,49] <- 0
    psi[d,t,62] <- 0
    psi[d,t,67] <- 0
    psi[d,t,80] <- 0
    } #d
    
    for (sf in 1:8){  
    for (s in 1:8){
    psi[(6*(sf-1))+1,t,(s*6-4)] <- phi[sf,1,t]*psi1[sf,s]
    psi[(6*(sf-1))+1,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+1,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+1,t,(s*6)] <- 0
    psi[(6*(sf-1))+2,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+2,t,(s*6-3)] <- phi[sf,2,t]*psi2[sf,s]
    psi[(6*(sf-1))+2,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+2,t,(s*6)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+3,t,(s*6-2)] <- phi[sf,3,t]*psi3[sf,s]
    psi[(6*(sf-1))+3,t,(s*6)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+4,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+4,t,(s*6)] <- phi[sf,4,t]*psi4[sf,s]
    psi[(6*(sf-1))+5,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+5,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+5,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+5,t,(s*6)] <- phi[sf,5,t]*psiAa[sf,s]
    psi[(6*(sf-1))+6,t,(s*6-4)] <- 0
    psi[(6*(sf-1))+6,t,(s*6-3)] <- 0
    psi[(6*(sf-1))+6,t,(s*6-2)] <- 0
    psi[(6*(sf-1))+6,t,(s*6)] <- phi[sf,6,t]*psiAb[sf,s]
    } #s
    
    psi[(6*(sf-1))+1,t,50] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*g9[t, 2]
    psi[(6*(sf-1))+1,t,51] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[(6*(sf-1))+1,t,52] <- phi[sf,1,t]*psi1[sf,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[(6*(sf-1))+1,t,53] <- phi[sf,1,t]*psi1[sf,9]*(1-p9[2,t])
    for (ex in c(54:61,63:66)){
    psi[(6*(sf-1))+1,t,ex] <- 0  }
    psi[(6*(sf-1))+1,t,68] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*g10[t, 2]
    psi[(6*(sf-1))+1,t,69] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[(6*(sf-1))+1,t,70] <- phi[sf,1,t]*psi1[sf,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[(6*(sf-1))+1,t,71] <- phi[sf,1,t]*psi1[sf,10]*(1-p10[2,t])
    for (ex in c(72:79,81:84)){
    psi[(6*(sf-1))+1,t,ex] <- 0  }
    
    for (ew in 50:53){
    psi[(6*(sf-1))+2,t,ew] <- 0}
    psi[(6*(sf-1))+2,t,54] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*g9[t, 3]
    psi[(6*(sf-1))+2,t,55] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[(6*(sf-1))+2,t,56] <- phi[sf,2,t]*psi2[sf,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[(6*(sf-1))+2,t,57] <- phi[sf,2,t]*psi2[sf,9]*(1-p9[3,t])
    for (er in c(58:61,63:66)){
    psi[(6*(sf-1))+2,t,er] <- 0  }
    for (ew in 68:71){
    psi[(6*(sf-1))+2,t,ew] <- 0}
    psi[(6*(sf-1))+2,t,72] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*g10[t, 3]
    psi[(6*(sf-1))+2,t,73] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[(6*(sf-1))+2,t,74] <- phi[sf,2,t]*psi2[sf,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[(6*(sf-1))+2,t,75] <- phi[sf,2,t]*psi2[sf,10]*(1-p10[3,t])
    for (er in c(76:79,81:84)){
    psi[(6*(sf-1))+2,t,er] <- 0  }
    
    for (ed in 50:57){
    psi[(6*(sf-1))+3,t,ed] <- 0}
    psi[(6*(sf-1))+3,t,58] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*g9[t, 4]
    psi[(6*(sf-1))+3,t,59] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[(6*(sf-1))+3,t,60] <- phi[sf,3,t]*psi3[sf,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[(6*(sf-1))+3,t,61] <- phi[sf,3,t]*psi3[sf,9]*(1-p9[4,t])
    for (ef in 63:66){
    psi[(6*(sf-1))+3,t,ef] <- 0 }
    for (ed in 68:75){
    psi[(6*(sf-1))+3,t,ed] <- 0}
    psi[(6*(sf-1))+3,t,76] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*g10[t, 4]
    psi[(6*(sf-1))+3,t,77] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[(6*(sf-1))+3,t,78] <- phi[sf,3,t]*psi3[sf,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[(6*(sf-1))+3,t,79] <- phi[sf,3,t]*psi3[sf,10]*(1-p10[4,t])
    for (ef in 81:84){
    psi[(6*(sf-1))+3,t,ef] <- 0 }
    
    for (es in 50:61){
    psi[(6*(sf-1))+4,t,es] <- 0}
    psi[(6*(sf-1))+4,t,63] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+4,t,64] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+4,t,65] <- phi[sf,4,t]*psi4[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+4,t,66] <- phi[sf,4,t]*psi4[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+4,t,es] <- 0}
    psi[(6*(sf-1))+4,t,81] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+4,t,82] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+4,t,83] <- phi[sf,4,t]*psi4[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+4,t,84] <- phi[sf,4,t]*psi4[sf,10]*(1-p10[6,t])
    
    for (es in 50:61){
    psi[(6*(sf-1))+5,t,es] <- 0}
    psi[(6*(sf-1))+5,t,63] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+5,t,64] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+5,t,65] <- phi[sf,5,t]*psiAa[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+5,t,66] <- phi[sf,5,t]*psiAa[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+5,t,es] <- 0}
    psi[(6*(sf-1))+5,t,81] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+5,t,82] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+5,t,83] <- phi[sf,5,t]*psiAa[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+5,t,84] <- phi[sf,5,t]*psiAa[sf,10]*(1-p10[6,t])
    
    for (es in 50:61){
    psi[(6*(sf-1))+6,t,es] <- 0}
    psi[(6*(sf-1))+6,t,63] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*g9[t, 6]
    psi[(6*(sf-1))+6,t,64] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+6,t,65] <- phi[sf,6,t]*psiAb[sf,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[(6*(sf-1))+6,t,66] <- phi[sf,6,t]*psiAb[sf,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[(6*(sf-1))+6,t,es] <- 0}
    psi[(6*(sf-1))+6,t,81] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*g10[t, 6]
    psi[(6*(sf-1))+6,t,82] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+6,t,83] <- phi[sf,6,t]*psiAb[sf,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[(6*(sf-1))+6,t,84] <- phi[sf,6,t]*psiAb[sf,10]*(1-p10[6,t])
    } #sf
    
    
    for (s in 1:8){
    psi[49,t,(s*6-4)] <- phi[9,1,t]*psi1[9,s]
    psi[49,t,(s*6-3)] <- 0
    psi[49,t,(s*6-2)] <- 0
    psi[49,t,(s*6)] <- 0
    for (i in 1:4){
    psi[49+i,t,(s*6-4)] <- 0
    psi[49+i,t,(s*6-3)] <- phi[9,2,t]*psi2[9,s]
    psi[49+i,t,(s*6-2)] <- 0
    psi[49+i,t,(s*6)] <- 0
    psi[53+i,t,(s*6-4)] <- 0
    psi[53+i,t,(s*6-3)] <- 0
    psi[53+i,t,(s*6-2)] <- phi[9,3,t]*psi3[9,s]
    psi[53+i,t,(s*6)] <- 0
    psi[57+i,t,(s*6-4)] <- 0
    psi[57+i,t,(s*6-3)] <- 0
    psi[57+i,t,(s*6-2)] <- 0
    psi[57+i,t,(s*6)] <- phi[9,4,t]*psi4[9,s] 
    psi[62+i,t,(s*6-4)] <- 0
    psi[62+i,t,(s*6-3)] <- 0
    psi[62+i,t,(s*6-2)] <- 0
    psi[62+i,t,(s*6)] <- phi[9,6,t]*psiAb[9,s]
    } #i
    psi[62,t,(s*6-4)] <- 0
    psi[62,t,(s*6-3)] <- 0
    psi[62,t,(s*6-2)] <- 0
    psi[62,t,(s*6)] <- phi[9,5,t]*psiAa[9,s]
    } #s
    
    for (s in 1:8){
    psi[67,t,(s*6-4)] <- phi[10,1,t]*psi1[10,s]
    psi[67,t,(s*6-3)] <- 0
    psi[67,t,(s*6-2)] <- 0
    psi[67,t,(s*6)] <- 0
    for (i in 1:4){
    psi[67+i,t,(s*6-4)] <- 0
    psi[67+i,t,(s*6-3)] <- phi[10,2,t]*psi2[10,s]
    psi[67+i,t,(s*6-2)] <- 0
    psi[67+i,t,(s*6)] <- 0
    psi[71+i,t,(s*6-4)] <- 0
    psi[71+i,t,(s*6-3)] <- 0
    psi[71+i,t,(s*6-2)] <- phi[10,3,t]*psi3[10,s]
    psi[71+i,t,(s*6)] <- 0
    psi[75+i,t,(s*6-4)] <- 0
    psi[75+i,t,(s*6-3)] <- 0
    psi[75+i,t,(s*6-2)] <- 0
    psi[75+i,t,(s*6)] <- phi[10,4,t]*psi4[10,s] 
    psi[80+i,t,(s*6-4)] <- 0
    psi[80+i,t,(s*6-3)] <- 0
    psi[80+i,t,(s*6-2)] <- 0
    psi[80+i,t,(s*6)] <- phi[10,6,t]*psiAb[10,s]
    } #i
    psi[80,t,(s*6-4)] <- 0
    psi[80,t,(s*6-3)] <- 0
    psi[80,t,(s*6-2)] <- 0
    psi[80,t,(s*6)] <- phi[10,5,t]*psiAa[10,s]
    } #s
    
    psi[49,t,50] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*g9[t, 2]
    psi[49,t,51] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[49,t,52] <- phi[9,1,t]*psi1[9,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[49,t,53] <- phi[9,1,t]*psi1[9,9]*(1-p9[2,t])
    for (ex in c(54:61,63:66)){
    psi[49,t,ex] <- 0  }
    
    for (i in 1:4){
    for (ew in 50:53){
    psi[49+i,t,ew] <- 0  }
    psi[49+i,t,54] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*g9[t, 3]
    psi[49+i,t,55] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[49+i,t,56] <- phi[9,2,t]*psi2[9,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[49+i,t,57] <- phi[9,2,t]*psi2[9,9]*(1-p9[3,t])
    for (er in c(58:61,63:66)){
    psi[49+i,t,er] <- 0  }
    for (ed in 50:57){
    psi[53+i,t,ed] <- 0}
    psi[53+i,t,58] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*g9[t, 4]
    psi[53+i,t,59] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[53+i,t,60] <- phi[9,3,t]*psi3[9,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[53+i,t,61] <- phi[9,3,t]*psi3[9,9]*(1-p9[4,t])
    for (ef in 63:66){
    psi[53+i,t,ef] <- 0 }
    for (es in 50:61){
    psi[57+i,t,es] <- 0}
    psi[57+i,t,63] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*g9[t, 6]
    psi[57+i,t,64] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[57+i,t,65] <- phi[9,4,t]*psi4[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[57+i,t,66] <- phi[9,4,t]*psi4[9,9]*(1-p9[6,t])
    for (es in 50:61){
    psi[62+i,t,es] <- 0}
    psi[62+i,t,63] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*g9[t, 6]
    psi[62+i,t,64] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62+i,t,65] <- phi[9,6,t]*psiAb[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62+i,t,66] <- phi[9,6,t]*psiAb[9,9]*(1-p9[6,t])
    } #i
    
    for (es in 50:61){
    psi[62,t,es] <- 0}
    psi[62,t,63] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*g9[t, 6]
    psi[62,t,64] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62,t,65] <- phi[9,5,t]*psiAa[9,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[62,t,66] <- phi[9,5,t]*psiAa[9,9]*(1-p9[6,t])
    
    
    psi[49,t,50+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*g10[t, 2]
    psi[49,t,51+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[49,t,52+18] <- phi[9,1,t]*psi1[9,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[49,t,53+18] <- phi[9,1,t]*psi1[9,10]*(1-p10[2,t])
    for (ex in c(54:61,63:66)){
    psi[49,t,ex+18] <- 0  }
    
    for (i in 1:4){
    for (ew in 50:53){
    psi[49+i,t,ew+18] <- 0  }
    psi[49+i,t,54+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*g10[t, 3]
    psi[49+i,t,55+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[49+i,t,56+18] <- phi[9,2,t]*psi2[9,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[49+i,t,57+18] <- phi[9,2,t]*psi2[9,10]*(1-p10[3,t])
    for (er in c(58:61,63:66)){
    psi[49+i,t,er+18] <- 0  }
    for (ed in 50:57){
    psi[53+i,t,ed+18] <- 0}
    psi[53+i,t,58+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*g10[t, 4]
    psi[53+i,t,59+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[53+i,t,60+18] <- phi[9,3,t]*psi3[9,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[53+i,t,61+18] <- phi[9,3,t]*psi3[9,10]*(1-p10[4,t])
    for (ef in 63:66){
    psi[53+i,t,ef+18] <- 0 }
    for (es in 50:61){
    psi[57+i,t,es+18] <- 0}
    psi[57+i,t,63+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*g10[t, 6]
    psi[57+i,t,64+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[57+i,t,65+18] <- phi[9,4,t]*psi4[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[57+i,t,66+18] <- phi[9,4,t]*psi4[9,10]*(1-p10[6,t])
    for (es in 50:61){
    psi[62+i,t,es+18] <- 0}
    psi[62+i,t,63+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*g10[t, 6]
    psi[62+i,t,64+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62+i,t,65+18] <- phi[9,6,t]*psiAb[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62+i,t,66+18] <- phi[9,6,t]*psiAb[9,10]*(1-p10[6,t])
    } #i
    
    for (es in 50:61){
    psi[62,t,es+18] <- 0}
    psi[62,t,63+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*g10[t, 6]
    psi[62,t,64+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62,t,65+18] <- phi[9,5,t]*psiAa[9,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[62,t,66+18] <- phi[9,5,t]*psiAa[9,10]*(1-p10[6,t])
    
    psi[67,t,50] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*g9[t, 2]
    psi[67,t,51] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[67,t,52] <- phi[10,1,t]*psi1[10,9]*p9[2,t]*g9[t, 2]*(1-g9[t, 2])/(1+g9[t, 2])
    psi[67,t,53] <- phi[10,1,t]*psi1[10,9]*(1-p9[2,t])
    for (ex in c(72:79,81:84)){
    psi[67,t,ex-18] <- 0  }
    
    for (i in 1:4){
    for (ew in 68:71){
    psi[67+i,t,ew-18] <- 0  }
    psi[67+i,t,72-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*g9[t, 3]
    psi[67+i,t,73-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[67+i,t,74-18] <- phi[10,2,t]*psi2[10,9]*p9[3,t]*g9[t, 3]*(1-g9[t, 3])/(1+g9[t, 3])
    psi[67+i,t,75-18] <- phi[10,2,t]*psi2[10,9]*(1-p9[3,t])
    for (er in c(76:79,81:84)){
    psi[67+i,t,er-18] <- 0  }
    for (ed in 68:75){
    psi[71+i,t,ed-18] <- 0}
    psi[71+i,t,76-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*g9[t, 4]
    psi[71+i,t,77-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[71+i,t,78-18] <- phi[10,3,t]*psi3[10,9]*p9[4,t]*g9[t, 4]*(1-g9[t, 4])/(1+g9[t, 4])
    psi[71+i,t,79-18] <- phi[10,3,t]*psi3[10,9]*(1-p9[4,t])
    for (ef in 81:84){
    psi[71+i,t,ef-18] <- 0 }
    for (es in 68:79){
    psi[75+i,t,es-18] <- 0}
    psi[75+i,t,81-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*g9[t, 6]
    psi[75+i,t,82-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[75+i,t,83-18] <- phi[10,4,t]*psi4[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[75+i,t,84-18] <- phi[10,4,t]*psi4[10,9]*(1-p9[6,t])
    for (es in 68:79){
    psi[80+i,t,es-18] <- 0}
    psi[80+i,t,81-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*g9[t, 6]
    psi[80+i,t,82-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80+i,t,83-18] <- phi[10,6,t]*psiAb[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80+i,t,84-18] <- phi[10,6,t]*psiAb[10,9]*(1-p9[6,t])
    } #i
    
    for (es in 68:79){
    psi[80,t,es-18] <- 0}
    psi[80,t,81-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*g9[t, 6]
    psi[80,t,82-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80,t,83-18] <- phi[10,5,t]*psiAa[10,9]*p9[6,t]*g9[t, 6]*(1-g9[t, 6])/(1+g9[t, 6])
    psi[80,t,84-18] <- phi[10,5,t]*psiAa[10,9]*(1-p9[6,t])
    
    
    psi[67,t,68] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*g10[t, 2]
    psi[67,t,69] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[67,t,70] <- phi[10,1,t]*psi1[10,10]*p10[2,t]*g10[t, 2]*(1-g10[t, 2])/(1+g10[t, 2])
    psi[67,t,71] <- phi[10,1,t]*psi1[10,10]*(1-p10[2,t])
    for (ex in c(72:79,81:84)){
    psi[67,t,ex] <- 0  }
    
    for (i in 1:4){
    for (ew in 68:71){
    psi[67+i,t,ew] <- 0  }
    psi[67+i,t,72] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*g10[t, 3]
    psi[67+i,t,73] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[67+i,t,74] <- phi[10,2,t]*psi2[10,10]*p10[3,t]*g10[t, 3]*(1-g10[t, 3])/(1+g10[t, 3])
    psi[67+i,t,75] <- phi[10,2,t]*psi2[10,10]*(1-p10[3,t])
    for (er in c(76:79,81:84)){
    psi[67+i,t,er] <- 0  }
    for (ed in 68:75){
    psi[71+i,t,ed] <- 0}
    psi[71+i,t,76] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*g10[t, 4]
    psi[71+i,t,77] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[71+i,t,78] <- phi[10,3,t]*psi3[10,10]*p10[4,t]*g10[t, 4]*(1-g10[t, 4])/(1+g10[t, 4])
    psi[71+i,t,79] <- phi[10,3,t]*psi3[10,10]*(1-p10[4,t])
    for (ef in 81:84){
    psi[71+i,t,ef] <- 0 }
    for (es in 68:79){
    psi[75+i,t,es] <- 0}
    psi[75+i,t,81] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*g10[t, 6]
    psi[75+i,t,82] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[75+i,t,83] <- phi[10,4,t]*psi4[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[75+i,t,84] <- phi[10,4,t]*psi4[10,10]*(1-p10[6,t])
    for (es in 68:79){
    psi[80+i,t,es] <- 0}
    psi[80+i,t,81] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*g10[t, 6]
    psi[80+i,t,82] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80+i,t,83] <- phi[10,6,t]*psiAb[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80+i,t,84] <- phi[10,6,t]*psiAb[10,10]*(1-p10[6,t])
    } #i
    for (es in 68:79){
    psi[80,t,es] <- 0}
    psi[80,t,81] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*g10[t, 6]
    psi[80,t,82] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80,t,83] <- phi[10,5,t]*psiAa[10,10]*p10[6,t]*g10[t, 6]*(1-g10[t, 6])/(1+g10[t, 6])
    psi[80,t,84] <- phi[10,5,t]*psiAa[10,10]*(1-p10[6,t])
    
    
    # Define probabilities of O(t)
    for (j in 1:8){
    for (cc in c(2:4,6)){
    po[((j-1)*6)+cc,t] <- p[j,cc,t] }
    po[1+(6*(j-1)),t] <- 0 #State that can only occur on first capture
    po[5+(6*(j-1)),t] <- 0 #State that can only occur on first capture
    }
    for (jj in c(50:52, 54:56, 58:60, 63:65, 68:70, 72:74, 76:78, 81:83)){
    po[jj,t] <- 1
    }
    po[49,t] <- 0 #State that can only occur on first capture
    po[62,t] <- 0 #State that can only occur on first capture
    po[53,t] <- 0 #Unobservable state
    po[57,t] <- 0 #Unobservable state
    po[61,t] <- 0 #Unobservable state
    po[66,t] <- 0 #Unobservable state
    po[67,t] <- 0 #State that can only occur on first capture
    po[80,t] <- 0 #State that can only occur on first capture
    po[71,t] <- 0 #Unobservable state
    po[75,t] <- 0 #Unobservable state
    po[79,t] <- 0 #Unobservable state
    po[84,t] <- 0 #Unobservable state
    
    # Calculate probability of non-encounter (dq) and reshape the array for the encounter probabilities      
    for (s in 1:ns){
    dp[s,t,s] <- po[s,t]
    dq[s,t,s] <- 1-po[s,t]
    } # s
    
    for (s in 1:(ns-1)){
    for (m in (s+1):ns){
    dp[s,t,m] <- 0
    dq[s,t,m] <- 0
    } # s
    } # m
    for (s in 2:ns){
    for (m in 1:(s-1)){
    dp[s,t,m] <- 0
    dq[s,t,m] <- 0
    } # s
    } # m
    } # t
    
    # Define the multinomial likelihood
    for (t in 1:((n.occasions-1)*ns)){
    marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
    }
    
    # Define the cell probabilities of the multistate m-array   
    # Define matrix U: product of probabilities of state-transition and non-encounter (this is just done because there is no product function for matrix multiplication in JAGS)
    for (t in 1:(n.occasions-2)){
    U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
    for (j in (t+1):(n.occasions-1)){
    U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,t,] %*% dq[,t,]
    }
    }
    U[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ones
    
    # Diagonal
    for (t in 1:(n.occasions-2)){
    pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*% dp[,t,]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% psi[,j,] %*% dp[,j,]
    }
    }
    pr[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- psi[,n.occasions-1,] %*% dp[,n.occasions-1,]
    
    # Below main diagonal
    for (t in 2:(n.occasions-1)){
    for (j in 1:(t-1)){
    pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
    } #j
    } #t
    
    # Last column: probability of non-recapture
    for (t in 1:((n.occasions-1)*ns)){
    pr[t,(n.occasions*ns-(ns-1))] <- 1-sum(pr[t,1:((n.occasions-1)*ns)])
    } #t
    }    ",fill = TRUE)
sink()

# Bundle data
n.unobs <- 8 #Number of unobservable states
n.notobs <- 15 #Number of observable states that dont appear anywhere in capture histories
ns <- length(unique(as.numeric(mat.enc3))) - 1 + n.unobs + n.notobs # calculate the number of states

jags.data <- list(marr = ms.arr, n.occasions = ncol(mat.enc3), rel = rowSums(ms.arr), ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), 
                  dista=distmat, ind=in., in.p=in.p, p.all.years=p.all.years)

inits <- function(){list(a.del=rep(0.5, length.out=3),   
                         mean.phi = 0.1, sigma.phi = runif(1, 0, 5),
                         mean.p = runif(2,min=0,max=1), sigma.p = runif(1, 0, 5), A.delp=rnorm(3),
                         p.both.enc9 = rnorm(4), p.both.enc10 = rnorm(1),
                         t.g9=rnorm(4), a.g.3 = 0.5,a.g.4 = 0.5,a.g.6 = 0.5, t.g10 = 0.5,
                         B=rnorm(n=1, mean=-0.3, sd=0.1),
                         theta2 = rnorm(1), theta3 = rnorm(1), theta4 = rnorm(1), thetaAa = rnorm(1), thetaAb = rnorm(1), repuls=rnorm(10), attract=rnorm(10)
)}  

# Parameters monitored
parameters <- c("phi.real", "sigma2.phi.real", "p.real", "sigma2.p.real", "p9.real","g9.real", "p10.real", "g10.real","B", "psi1", "psi2","psi3", "psiAa", "psiAb")

# MCMC settings
ni <- 15 #00 
nt <- 1
nb <- 2 #00
nc <- 3

# Call JAGS from R, using jagsUI syntax
starttime <- Sys.time()
m11 <- jags(jags.data, inits, parameters, "model11.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=FALSE)
endtime <- Sys.time()
options(max.print=999999) #Default is 99999
print(m11, digits = 3) # 
print(m11$summary[1:400,], digits = 3)
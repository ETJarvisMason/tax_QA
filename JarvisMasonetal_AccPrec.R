# Jarvis Mason, ET
# Quantify accuracy and precision by species, across taxonomists
# and separately, by species and taxonomist
# 2023


# load libraries ----------------------------------------------------------


library(R2jags)
library(MCMCvis)
library(dplyr)
library(tidyr)

setwd('C:/Users/etmason/Documents/CASeaGrant2018/Morphology/ASawkins')

# you may wish to use more than one training data set (here, we use two batches of data)
# load training data (Batch 1) ---------------------------------------------------------------

# read in training data set
batch1 <- read.csv("batch1_seq.csv")

# put it into long format
# here, morphID refers to the taxonomist id
df <- gather(batch1,key="taxonomist",value = "morphID",-"molecular.ID")
names(df)

# batch 1 data
b1 <- df %>% 
  dplyr::rename(molecID = "molecular.ID") %>% # known species identification
  mutate(success = 
           case_when(
             morphID == molecID ~ 1,
             morphID != molecID ~ 0) # select cases where the ID was correct
           )
b1$taxonomist <- as.numeric(as.factor(b1$taxonomist)) # change names to numbers 
       

# Format data for accuracy ----------------------------------------------------------------

# summarize data by species, across all taxonomists 

bsb.b1 <- b1 %>% 
  filter(molecID == "BSB") %>% 
  group_by(molecID) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())
kb.b1 <- b1 %>% 
  filter(molecID == "KB") %>% 
  group_by(molecID) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())


# summarize data by species and taxonomist
bsb.b1.tax <- b1 %>% 
  filter(molecID == "BSB") %>% 
  group_by(molecID, taxonomist) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

kb.b1.tax <- b1 %>% 
  filter(molecID == "KB") %>% 
  group_by(molecID, taxonomist) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

# Format data for precision ----------------------------------------------------------------

# summarize data by species, across all taxonomists 

bsb.b1.prec <- b1 %>% 
  filter(morphID == "BSB") %>% 
  group_by(morphID) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

kb.b1.prec <- b1 %>% 
  filter(morphID == "KB") %>% 
  group_by(morphID) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

# summarize data by species and taxonomist
bsb.b1.prectax <- b1 %>% 
  filter(morphID == "BSB") %>% 
  group_by(morphID,taxonomist) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

kb.b1.prectax <- b1 %>% 
  filter(morphID == "KB") %>% 
  group_by(morphID, taxonomist) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

# models to estimate probability of success
# we are going to use this two ways:
#1) to model accuracy (What is the probability that any species identification is correct)
#2) to model precision (What is the probability that a specific species positive identification (morph ID) is a true positive?
#Precision = TP / (TP + FP)

# Accuracy and precision of species identifications across all taxonomists
# Species Model ------------------------------------------------------

# Tests for Accuracy --------------------------------------------------------------

# we will use this model to generate estimates of theta by species 

sink("species_mod.jags")
cat("
    model{
    
    # Likelihood
    Y ~ dbinom(theta,n)
    
    # Prior
    theta ~ dbeta(a, b)
    }
    ",fill = TRUE)
sink()

parameters <- c("theta") 

# model specifications regarding how many iterations, etc.
ni <- 40000
nt <- 40
nb <- 20000
nc <- 3


set.seed(123) 

####BSB######
Y <- bsb.b1$TruePositive
n <- bsb.b1$Samples
# this defines the data
jags.data <- list(Y = Y,n = n,a = a,b = b)
b1.BSB <- jags(jags.data, inits=NULL, parameters, 
               "species_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
# print results
print(b1.BSB, dig = 3)

####KB######
# Data for Kelp Bass - overall
Y <- kb.b1$TruePositive
n <- kb.b1$Samples
# this defines the data
jags.data <- list(Y = Y,
                  n = n,
                  a = a,
                  b = b)
b1.KB <- jags(jags.data, inits=NULL, parameters, 
              "species_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd()) 
print(b1.KB, dig = 3)
# plot trace plots and posterior distributions

MCMCplot(b1.BSB,
         b1.KB,
         params = 'theta',
         ref = 0.5,
         col = "dark blue",
         col2 = "dark orange",
         main = 'Batch 1 Classification Accuracy',
         xlim = c(0.4,1.05),
         xlab = "Probability",
         horiz = TRUE,
         labels = NULL,
         offset = 0.3)

# Tests for Precision ------------------------------------------------------


######BSB######
Y <- bsb.b1.prec$TruePositive
n <- bsb.b1.prec$Samples
jags.data <- list(Y = Y,n = n,a = a,b = b)
b1.BSBP <- jags(jags.data, inits=NULL, parameters, 
                "species_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(b1.BSBP, dig = 3)

########KB#######
Y <- kb.b1.prec$TruePositive
n <- kb.b1.prec$Samples
jags.data <- list(Y = Y,n = n,a = a,b = b)
b1.KBP <- jags(jags.data, inits=NULL, parameters, 
               "species_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(b1.KBP, dig = 3)

MCMCplot(b1.BSBP,
         b1.KBP,
         params = 'theta',
         ref = 0.5,
         col = "dark blue",
         col2 = "dark orange",
         main = 'Batch 1 Classification Precision',
         xlim = c(0.4,1.05),
         xlab = "Probability",
         labels = NULL,
         offset = 0.3)



# Taxonomist Model (species by taxonomist) --------------------------------

# we will use this model to generate estimates of theta by species and taxonomist 
# here, taxonomist is treated as a random effect
# if there are fewer than five taxonomists, then taxonomist
# needs to be treated as a fixed effect...see fixed effect model below

sink("spptaxon_mod.jags")
cat("
    model{
    
    # Likelihood
    
    for (i in 1:length(Y)){ # for every taxonomist do the following
    
    x[i]~dnorm(mu,tau)

    theta[i] <- 1/(1 + exp(-x[i]))
    
    Y[i] ~ dbinom(theta[i], n[i])
    
    }# end for loop


    # Prior
    
    mu ~ dnorm(0, 0.001) # hyperprior
    

    sigma ~ dunif(0, 10) # put a uniform on the sd

    tau <- pow(sigma, -2) # precision (specific to JAGS)
    
    # for (i in 1:length(Y)){ # for every taxonomist do the following
    # theta[i] ~ dbeta(a, b)
    # }
    }
    ",fill = TRUE)
sink()

parameters <- c("mu","theta","sigma") 

# model specifications regarding how many iterations, etc.
ni <- 40000
nt <- 40
nb <- 20000
nc <- 3


# Tests for Accuracy ------------------------------------------------------


####BSB######
Y <- bsb.b1.tax$TruePositive
n <- bsb.b1.tax$Samples
# jags.data <- list(Y = Y,n = n,a = a,b = b) # fixed effects
jags.data <- list(Y = Y,n = n)# random effects
b1.BSB.tax <- jags(jags.data, inits=NULL, parameters, 
                   "spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(b1.BSB.tax, dig = 3)
####KB######
Y <- kb.b1.tax$TruePositive
n <- kb.b1.tax$Samples
jags.data <- list(Y = Y,n = n)
b1.KB.tax <- jags(jags.data, inits=NULL, parameters, 
                  "spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(b1.KB.tax, dig = 3)
b1.acc <- MCMCplot(b1.BSB.tax,
         b1.KB.tax,
         params = 'theta',
         ref = 0.5,
         col = "dark blue",
         col2 = "dark orange",
         main = 'Batch 1 Taxonomist Accuracy',
         xlim = c(0.4,1.05),
         xlab = "Probability",
         labels = c('1','2','3','4','5'),
         offset = 0.1)

# Tests for Precision ------------------------------------------------------

#######BSB#####
Y <- bsb.b1.prectax$TruePositive
n <- bsb.b1.prectax$Samples
jags.data <- list(Y = Y,n = n)
b1.BSB.tax.prec <- jags(jags.data, inits=NULL, parameters, 
                   "spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(b1.BSB.tax.prec, dig = 3)
#######KB########
Y <- kb.b1.prectax$TruePositive
n <- kb.b1.prectax$Samples
jags.data <- list(Y = Y,n = n,a = a,b = b)
b1.KB.tax.prec <- jags(jags.data, inits=NULL, parameters, 
                  "spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(b1.KB.tax.prec, dig = 3)
b1.prec <-MCMCplot(b1.BSB.tax.prec,
         b1.KB.tax.prec,
         params = 'theta',
         ref = 0.5,
         col = "dark blue",
         col2 = "dark orange",
         main = 'Batch 1 Taxonomist Precision',
         xlim = c(0.4,1.05),
         xlab = "Probability",
         labels = c('1','2','3','4','5'),
         offset = 0.1)


# load training data (Batch 2) ---------------------------------------------------------------

# read in training data set
batch2 <- read.csv("Batch2_final_sequences.csv", stringsAsFactors = FALSE)
df <- batch2
names(df)
b2 <- df %>% 
  dplyr::rename(molecID = "Molecular.ID",
                taxonomist = "Sorter.Initials",
                morphID = "Sorter.ID") %>%
  mutate(success = 
           case_when(
             morphID == molecID ~ 1,
             morphID != molecID ~ 0) # select cases where the ID was correct
  )
b2$taxonomist <- as.numeric(as.factor(b2$taxonomist)) # change names to numbers 
b2  

bsb.b2 <- b2 %>% 
  filter(molecID == "BSB") %>% 
  group_by(molecID) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())
kb.b2 <- b2 %>% 
  filter(molecID == "KB") %>% 
  group_by(molecID) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

# Format data for accuracy ----------------------------------------------------------------

# summarize data by species, across all taxonomists 
bsb.b2.tax <- b2 %>% 
  filter(molecID == "BSB") %>% 
  group_by(molecID, taxonomist) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

# summarise data by species and taxonomist
kb.b2.tax <- b2 %>% 
  filter(molecID == "KB") %>% 
  group_by(molecID, taxonomist) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

# Format data for precision ----------------------------------------------------------------

# summarize data by species, across all taxonomists 

bsb.b2.prec <- b2 %>% 
  filter(morphID == "BSB") %>% 
  group_by(morphID) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

kb.b2.prec <- b2 %>% 
  filter(morphID == "KB") %>% 
  group_by(morphID) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

# summarize data by species and taxonomist 

bsb.b2.prectax <- b2 %>% 
  filter(morphID == "BSB") %>% 
  group_by(morphID,taxonomist) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

kb.b2.prectax <- b2 %>% 
  filter(morphID == "KB") %>% 
  group_by(morphID, taxonomist) %>% 
  dplyr::summarize(TruePositive =sum(success),Samples = n())

# Accuracy and precision of species identifications across all taxonomists
# Species Model ------------------------------------------------------

# Tests for Accuracy --------------------------------------------------------------

# we will use this model to generate estimates of theta by species 

sink("species_mod.jags")
cat("
    model{
    
    # Likelihood
    Y ~ dbinom(theta,n)
    
    # Prior
    theta ~ dbeta(a, b)
    }
    ",fill = TRUE)
sink()


parameters <- c("theta") 

# model specifications regarding how many iterations, etc.
ni <- 40000
nt <- 40
nb <- 20000
nc <- 3

set.seed(123) 

####BSB######
Y <- bsb.b2$TruePositive
n <- bsb.b2$Samples
# this defines the data
jags.data <- list(Y = Y,n = n,a = a,b = b)
b2.BSB <- jags(jags.data, inits=NULL, parameters, 
               "species_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
# print results
print(b2.BSB, dig = 3)

####KB######
# Data for Kelp Bass - overall
Y <- kb.b2$TruePositive
n <- kb.b2$Samples
# this defines the data
jags.data <- list(Y = Y,
                  n = n,
                  a = a,
                  b = b)
b2.KB <- jags(jags.data, inits=NULL, parameters, 
              "species_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd()) 
print(b2.KB, dig = 3)
# plot trace plots and posterior distributions
MCMCplot(b2.BSB,
         b2.KB,
         params = 'theta',
         ref = 0.5,
         col = "dark blue",
         col2 = "dark orange",
         main = 'Batch 2 Classification Accuracy',
         xlim = c(0,1.05),
         xlab = "Probability",
         offset = 0.3)

# Tests for Precision ------------------------------------------------------

######BSB######
Y <- bsb.b2.prec$TruePositive
n <- bsb.b2.prec$Samples
jags.data <- list(Y = Y,n = n,a = a,b = b)
b2.BSBP <- jags(jags.data, inits=NULL, parameters, 
                "species_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(b2.BSBP, dig = 3)

######KB######
Y <- kb.b2.prec$TruePositive
n <- kb.b2.prec$Samples
jags.data <- list(Y = Y,n = n,a = a,b = b)
b2.KBP <- jags(jags.data, inits=NULL, parameters, 
               "species_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(b2.KBP, dig = 3)
MCMCplot(b2.BSBP,
         b2.KBP,
         params = 'theta',
         ref = 0.5,
         col = "dark blue",
         col2 = "dark orange",
         main = 'Batch 2 Classification Precision',
         xlim = c(0,1.05),
         xlab = "Probability",
         offset = 0.3)

# Taxonomist Model (species by taxonomist) --------------------------------
# Code below was used for Manuscript Figure 2

# we will use this model to generate estimates of accuracy and precision by species and taxonomist 
# here, taxonomist is treated as a fixed effect

####Fixed Effects Model##########
sink("spptaxon_mod.jags")
cat("
    model{
    
    # Likelihood
    for (i in 1:length(Y)){ # for every taxonomist do the following
    
    Y[i] ~ dbinom(theta[i],n[i])
    }
    # Prior
    for (i in 1:length(Y)){ # for every taxonomist do the following
    theta[i] ~ dbeta(a, b)
    }
    }
    ",fill = TRUE)
sink()


# here, taxonomist is treated as a random effect
####Random Effects Model####
sink("re_spptaxon_mod.jags")
cat("
    model{
    
    # Likelihood
    
    for (i in 1:length(Y)){ # for every taxonomist do the following

    # x[i]~dnorm(mu,tau)
    x[i] ~ dnorm(logit.mu,tau)
    theta[i] <- 1/(1 + exp(-x[i]))
    Y[i] ~ dbinom(theta[i], n[i])
    }# end for loop

    # Prior
    
    # mu ~ dnorm(0, 0.001) # hyperprior PERILOUS...result ends up being u-shaped
    mu ~ dbeta(1,1) # global mean
    logit.mu <- log(mu)/(1 - log(mu))
    
    sigma ~ dunif(0, 10) # put a uniform on the sd
    tau <- pow(sigma, -2) # precision (specific to JAGS)
    
    # for (i in 1:length(Y)){ # for every taxonomist do the following
    # theta[i] ~ dbeta(a, b)
    # }
    }
    ",fill = TRUE)
sink()


parameters <- c("theta") #fixed effects model
parameters <- c("mu", "theta", "sigma") #random effects

# model specifications regarding how many iterations, etc.
ni <- 40000
nt <- 40
nb <- 20000
nc <- 3


# Tests for Accuracy --------------------------------------------------------------

# we will use this model to generate estimates of theta by species and taxonomist 

####BSB######
Y <- bsb.b2.tax$TruePositive
n <- bsb.b2.tax$Samples

# use for Fixed Effects Model
jags.data <- list(Y = Y,n = n,a = a,b = b) #for Fixed Effects Model
b2.BSB.tax <- jags(jags.data, inits=NULL, parameters, 
                   "spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# use for Random Effects Model
jags.data <- list(Y = Y,n = n) #for Random Effects Model

b2.BSB.tax <- jags(jags.data, inits=NULL, parameters, 
                   "re_spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(b2.BSB.tax, dig = 3)
####KB######
Y <- kb.b2.tax$TruePositive
n <- kb.b2.tax$Samples

# use for Fixed Effects Model
jags.data <- list(Y = Y,n = n,a = a,b = b) #for Fixed Effects Model
b2.KB.tax <- jags(jags.data, inits=NULL, parameters, 
                  "spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# use for Random Effects Model
jags.data <- list(Y = Y,n = n) #for Random Effects Model
b2.KB.tax <- jags(jags.data, inits=NULL, parameters, 
                  "re_spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


print(b2.KB.tax, dig = 3)

b2.acc <- MCMCplot(b2.BSB.tax,
                   b2.KB.tax,
                   params = 'theta',
                   ref = 0.5,
                   col = "dark blue",
                   col2 = "dark orange",
                   main = 'Batch 2 Taxonomist Accuracy',
                   xlim = c(0,1.05),
                   xlab = "Probability",
                   labels = c('1','2','3'),
                   offset = 0.1)
# Tests for Precision --------------------------------------------------------------

# we will use this model to generate estimates of theta by species and taxonomist 

######BSB######
Y <- bsb.b2.prectax$TruePositive
n <- bsb.b2.prectax$Samples

# use for Fixed Effects Model
jags.data <- list(Y = Y,n = n,a = a,b = b) #for Fixed Effects Model
b2.BSB.tax.prec <- jags(jags.data, inits=NULL, parameters, 
                        "spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
# use for Random Effects Model
jags.data <- list(Y = Y,n = n) #for Random Effects Model
b2.BSB.tax.prec <- jags(jags.data, inits=NULL, parameters, 
                        "re_spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(b2.BSB.tax.prec, dig = 3)

######BSB######
Y <- kb.b2.prectax$TruePositive
n <- kb.b2.prectax$Samples

# use for Fixed Effects Model
jags.data <- list(Y = Y,n = n,a = a,b = b) #for Fixed Effects Model
b2.KB.tax.prec <- jags(jags.data, inits=NULL, parameters, 
                       "spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# use for Random Effects Model
jags.data <- list(Y = Y,n = n) #for Random Effects Model
b2.KB.tax.prec <- jags(jags.data, inits=NULL, parameters, 
                       "re_spptaxon_mod.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(b2.KB.tax.prec, dig = 3)
b2.prec <- MCMCplot(b2.BSB.tax.prec,
                    b2.KB.tax.prec,
                    params = 'theta',
                    ref = 0.5,
                    col = "dark blue",
                    col2 = "dark orange",
                    main = 'Batch 2 Taxonomist Precision',
                    xlim = c(0,1.05),
                    xlab = "Probability",
                    labels = c('1','2','3'),
                    offset = 0.1)

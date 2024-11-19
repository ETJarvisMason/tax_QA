# 
# Explore morphological characters for species delineation and
# identify the most important characters
# 2023


# load libraries ----------------------------------------------------------


library(tidyverse)
library(readr)
library(R2jags)
library(mcmcplots)
library(stringr)
library(haven) # remove attributes
library(ggmcmc)
library(ggtext)
library(MCMCvis)
library(patchwork)
library(randomForest)
library(party)


theme_Publication <- function(base_size=20, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2),vjust = 0), #hjust = 0.5),
            # text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border =  element_blank(),#element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),#element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",# bottom
            legend.direction = "vertical",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(5,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}



# load data ---------------------------------------------------------------

# file with molecular identities
id <- read.csv(file = "ids.csv", stringsAsFactors = FALSE)
id <- id[, c(2, 4)] # get rid of extraneous columns

# remove specimens with unknown identities
id <- id %>%
  drop_na() %>%
  rename(Label = label)

# taxonomist blind morphological characters table

one <- read.csv(file = "morphtab_1.csv")
two <- read.csv(file = "morphtab_2.csv")


# explore utility of characters for species delineation ------------

# define function to run multinomial logistic regression model
# this function requires the user to specify the data (sorter/taxonomist)
# and any separate category of characters (stage, sex, etc.)
# for the purpose of this example, we include only two stages
run_anal <- function(sorter, stage) {
  # assign taxonomist blind morphological characters table
  # alternatively, the tables can be merged and a column for sorter added
  # for filtering
  
  if (sorter == "1") {
    tab <- one
  } else {
    tab <- two
  }
  tab
  # join characters table with known species ids
  tab <- left_join(tab, id, "Label")
  
  # here, "stage" refers to larval developmental stage, but this column
  # could theoretically be renamed to "category" to reflect any category of interest
  # (life stage, sex, region, etc.)
  # if there is no reason to split the characters by category, then the column "stage"
  # can be filled with a single generic number or letter
  
  
  if (stage == "PrF") {
    # here, we filter the characters table by stage (category)
    # preflexion larvae
    prf <- tab
    prf <- tab %>%
      filter(stage == "PrF") %>%
      mutate_if(sapply(prf, is.character), as.factor) %>%
      select(vpatch, numv, locv, lightv, dnumloc, species) %>%  # input stage-specific characters for considertion here
      drop_na()
    
    prf <- as.data.frame(prf)
    
    prf <- prf[, colSums(is.na(prf)) == 0]
    
    numspp <- length(unique(prf$species))

  #-------model in JAGS --------------
    species.mod <- function()  {
      for (i in 1:N) {
        species[i] ~ dcat(p[i, 1:J])
        
        for (j in 1:J) {
          log(q[i, j]) <-  b[1, j] +
            b[2, j] * vpatch[i] +
            b[3, j] * numv[i] +
            b[4, j] * locv[i] +
            b[5, j] * lightv[i] +
            b[6, j] * dnumloc[i]
          # s
          p[i, j] <-
            q[i, j] / sum(q[i, 1:J])  ## should be familiar from MLE notes: q is exp(Xb)
        }   # close J loop
      }  # close N loop
      
      # how many beta parameters? Equals number of characters + 1
      for (k in 1:6) {
        b[k, 1] <-
          0 ## MUST set the first set of covariates (for the first outcome category) to 0
        for (j in 2:J) {
          b[k, j] ~ dnorm(0, 0.1)
        }  # close J loop
      }  # close K loop
    }  # close model loop
    
    species.params <- c("b", "p")
    
    species.inits <- function() {
      list(b = matrix(
        c(rep(NA, 6),
          rep(0, 6*(numspp-1))),
        nrow = 6,
        ncol = numspp,
        byrow = FALSE
      ))
    }
    
    dprf <-
      as.list(prf)
    dprf$N <- nrow(prf)
    dprf$J <- length(as.numeric(levels(as.factor(prf$species))))
    
    species.fit <-
      jags(
        data = dprf,
        inits = species.inits,
        species.params,
        n.chains = 3,
        n.iter = 40000,
        n.burnin = 11000,
        n.thin = 10,
        model.file = species.mod
      )
    print(species.fit, dig = 3)
    
    tbl <- prf
  }#PreFlexion
  
  
  if (stage == "F") {
    # Flexion larvae
    f <- tab
    f <- tab %>%
      dplyr::filter(stage == "F") %>%
      mutate_if(sapply(f, is.character), as.factor) %>%
      select(vpatchF, vchar, crown, pec, dnumlocF, species) %>% # select for characters of interest and include the known identity column
      drop_na()
    
    numspp <- length(unique(f$species))
    
    #-------model in JAGS --------------
    species.mod <- function()  {
      for (i in 1:N) {
        species[i] ~ dcat(p[i, 1:J])
        
        for (j in 1:J) {
          log(q[i, j]) <-  b[1, j] +
            b[2, j] * vpatchF[i] +
            b[3, j] * vchar[i] +
            b[4, j] * crown[i] +
            b[5, j] * pec[i] +
            b[6, j] * dnumlocF[i]
          # s
          p[i, j] <-
            q[i, j] / sum(q[i, 1:J])  ## should be familiar from MLE notes: q is exp(Xb)
        }   # close J loop
      }  # close N loop
      
    # how many beta parameters? Equals number of characters + 1
        for (k in 1:6) {
        b[k, 1] <-
          0 ## MUST set the first set of covariates (for the first outcome category) to 0
        for (j in 2:J) {
          b[k, j] ~ dnorm(0, 0.1)
        }  # close J loop
      }  # close K loop
    }  # close model loop
    
    species.params <- c("b", "p")
    
     species.inits <- function() {
      list(b = matrix(
        c(rep(NA, 6),
          rep(0, 6*(numspp-1))),
        nrow = 6,
        ncol = numspp,
        byrow = FALSE
      ))
    }
    
    df <- as.list(f)
    df$N <- nrow(f)
    df$J <- length(as.numeric(levels(as.factor(f$species))))
    
    species.fit <-
      jags(
        data = df,
        inits = species.inits,
        species.params,
        n.chains = 3,
        n.iter = 40000,
        n.burnin = 11000,
        n.thin = 10,
        model.file = species.mod
      )
    print(species.fit, dig = 3)
    
    tbl <- f
  } # Flexion
  


# Plots

# for manuscript Figure 4  
mod <- species.fit

# Labels for coefficient plot (needs to be modified depending on number
# of traits/characters and species)  

# labels = c(
#     'Trait 1, Sp A',
#     'Trait 2, Sp A',
#     'Trait 3, Sp A',
#     'Trait 4, Sp A',
#     'Trait 5, Sp A',
#     'Trait 6, Sp A',
#     'Trait 1, Sp B',
#     'Trait 2, Sp B',
#     'Trait 3, Sp B',
#     'Trait 4, Sp B',
#     'Trait 5, Sp B',
#     'Trait 6, Sp B',
#     'Trait 1, Sp C',
#     'Trait 2, Sp C',
#     'Trait 3, Sp C',
#     'Trait 4, Sp C',
#     'Trait 5, Sp C',
#     'Trait 6, Sp C'
#   )
  
  MCMCplot(
    mod,
    params = 'b',
    ci = c(50, 90),
    guide_lines = TRUE,
    HPD = TRUE,
    ref = NULL,
    ref_ovl = TRUE
  )
  
  # for manuscript Figure 3   
  # define data for ggplot
  d.vec <- mod[["BUGSoutput"]][["summary"]]
  
  d <-
    as.data.frame(d.vec[which(rownames(d.vec) == "p[1,1]"):nrow(d.vec),])
  colnames(d) <-
    c("mean",
      "sd",
      "lCI",
      "25.CI",
      "50.CI",
      "75.CI",
      "uCI",
      "Rhat",
      "n.eff")
  d$num <- seq(1, nrow(tbl) * numspp, 1)
  
  d$species <- NA
  d$species[1:nrow(tbl)] <- "P. clathratus"
  
  if (numspp == 3) {
    d$species[nrow(tbl) + 1:nrow(tbl)] <- "P. maculatofasciatus"
    d$species[(nrow(tbl) * 2) + 1:nrow(tbl)] <- "P. nebulifer"
    d$ind <- rep(1:nrow(tbl), numspp)
  }  else{
    d$species[nrow(tbl) + 1:nrow(tbl)] <- "P. nebulifer"
    d$ind <- rep(1:nrow(tbl), numspp)
  }
  
  # add true IDs
  
  nrow <- nrow(d)
  d[nrow + 1:nrow(tbl),] <- NA
  d$mean[nrow + 1:nrow(tbl)] <- 1.05
  d$species[nrow + 1:nrow(tbl)] <- as.character(tbl$species)
  d[nrow + 1:nrow(tbl), c(3, 7)] <- 1.05
  d$ind <- rep(1:nrow(tbl), numspp + 1)
  
  if (sorter == "1") {
    plot <-   ggplot() +
      geom_point(data = d,
                 aes(x = ind, y = mean, color = species),
                 group = "species") +
      geom_errorbar(data = d, aes(
        x = ind,
        ymin = lCI,
        ymax = uCI,
        color = species
      )) +
      {
        if (numspp == 3)
          scale_color_manual(
            values = c("#00AFBB", "#E7B800", "#FC4E07"),
            labels = c(
              '*P. clathratus*',
              '*P. maculatofasciatus*',
              '*P. nebulifer*'
            )
          )
        else
          scale_color_manual(
            values = c("#00AFBB", "#FC4E07"),
            labels = c('*P. clathratus*', '*P. nebulifer*')
          )
      } +
      scale_y_continuous(
        name = "Probability",
        limits = c(0, 1.05),
        breaks = seq(0, 1.05, .1)
      ) +
      labs(x = "Individual larvae") +
      theme_Publication() +
      theme(legend.position = "none")
  } else{
    plot <-   ggplot() +
      geom_point(data = d,
                 aes(x = ind, y = mean, color = species),
                 group = "species") +
      geom_errorbar(data = d, aes(
        x = ind,
        ymin = lCI,
        ymax = uCI,
        color = species
      )) +
      {
        if (numspp == 3)
          scale_color_manual(
            values = c("#00AFBB", "#E7B800", "#FC4E07"),
            labels = c(
              '*P. clathratus*',
              '*P. maculatofasciatus*',
              '*P. nebulifer*'
            )
          )
        else
          scale_color_manual(
            values = c("#00AFBB", "#FC4E07"),
            labels = c('*P. clathratus*', '*P. nebulifer*')
          )
      } +
      scale_y_continuous(
        name = "Probability",
        limits = c(0, 1.05),
        breaks = seq(0, 1.05, .1)
      ) +
      labs(x = "Individual larvae") +
      theme_Publication() +
      theme(legend.text = element_markdown(),
            legend.title = element_blank())
  }
  
  plot
}#end function



# plot results of multinomial logistic regression model ------------------------------------------------------



PrF_results_1 <- run_anal(sorter = "1", stage = "PrF")
PrF_results_2 <- run_anal(sorter = "2", stage = "PrF")

F_results_1 <- run_anal(sorter = "1", stage = "F")
F_results_2 <- run_anal(sorter = "2", stage = "F")

# print plots together
a <- PrF_results_1[["plot_env"]][["plot"]]
b <- PrF_results_2[["plot_env"]][["plot"]]

c <- F_results_1[["plot_env"]][["plot"]]
d <- F_results_2[["plot_env"]][["plot"]]

# Manuscript Figure 3 (model results for both taxonomists)
a / b + plot_annotation(
  tag_levels = 'A',
  theme = theme(
    plot.title = element_text(size = 24),
    axis.title = element_text(face = "bold", size = rel(0.5))
  )
)



# identity most important characters --------------------------------------


# set random seed to make results reproducible:
set.seed(17)
# insert clean data here - this is just data for one taxonomist
rf <- as.data.frame(PrF_results_1[["plot_env"]][["prf"]])   # select data from above
dat <- rf 

# Generate a random sample of different lengths to create training and validation sets
index1 <- sample(1:nrow(dat), size = nrow(dat))
index2 <- sample(1:nrow(dat), size = nrow(dat) / 1.25)

# Assign the species identities to the training and validation data sets
training <- dat [index1, 1:ncol(dat)]
training$species <- factor(training$species)
validation <- dat [index2, 1:ncol(dat)]
validation$species <- factor(validation$species)

# Perform training:
rf_classifier = randomForest(
  species ~ .,
  data = training,
  ntree = 500,
  mtry = 2,
  importance = TRUE
)

rf_classifier

# plots results of variable importance ------------------------------------



varImpPlot(rf_classifier)

# fancier plot
create_rfplot <- function(rf, type) {
  imp <- importance(rf, type = type, scale = F)
  
  featureImportance <-
    data.frame(Feature = row.names(imp), Importance = imp[, 1])
  
  p <-
    ggplot(featureImportance, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity",
             fill = "#53cfff",
             width = 0.65) +
    coord_flip() +
    theme_light(base_size = 20) +
    theme(
      axis.title.x = element_text(size = 15, color = "black"),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = 15, color = "black"),
      axis.text.y  = element_text(size = 15, color = "black")
    )
  return(p)
}

create_rfplot(rf_classifier, type = 1)

# additional ways of assessing the overall performance of the characters
# validation set assessment #1: looking at confusion matrix
prediction_for_table <-
  predict(rf_classifier, validation[,-ncol(validation)])

observed <- as.data.frame(validation)

table(observed = observed[, ncol(observed)], predicted = prediction_for_table)

# Validation set assessment #2: ROC curves and AUC

# Need to import ROCR package for ROC curve plotting:
library(ROCR)

# Calculate the probability of new observations belonging to each class

prediction_for_roc_curve <-
  predict(rf_classifier, validation[,-ncol(validation)], type = "prob")

# Use pretty colours:
pretty_colours <- c("#F8766D", "#00BA38", "purple")
# Specify the different classes
classes <- levels(validation$species)
# For each class
for (i in 1:length(classes)) {
  # Define which observations belong to class[i]
  true_values <-
    ifelse(observed[, ncol(observed)] == classes[i], 1, 0)
  # Assess the performance of classifier for class[i]
  pred <- prediction(prediction_for_roc_curve[, i], true_values)
  perf <- performance(pred, "tpr", "fpr")
  if (i == 1)
  {
    plot(perf, main = "ROC Curve", col = pretty_colours[i])
  }
  else
  {
    plot(perf,
         main = "ROC Curve",
         col = pretty_colours[i],
         add = TRUE)
  }
  # Calculate the AUC and print it to screen
  auc.perf <- performance(pred, measure = "auc")
  
  print(auc.perf@y.values)
}


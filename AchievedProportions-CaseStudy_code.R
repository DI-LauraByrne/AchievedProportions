##Install Packages #############################################################################################################################################

library(DImodels) # For creating interaction columns
library(nlme) # For gls regression
library(ggplot2) # For stacked barcharts
library(Compositional) # ZADR
library(dplyr) # Data manipulation

################################################################################################################################################################


##CONTINUOUS ACHIEVED DI MODEL
##Data Preparation #############################################################################################################################################

setwd("C:/Users/Administrator/OneDrive - TCDUD.onmicrosoft.com/Documents/PhD Stuff/Tasks/Task 11 - Achieved Proportion Modelling/Cedar Creek Data")

BioCON <- read.csv("e141_Plant aboveground biomass data.csv") # biomass response data
SetUp <- read.csv("e141_PlotSown.csv") # experiment design (species compositions)

SetUp <- SetUp[-372, -19:-22] # remove unneeded columns

SetUp <- SetUp[-which(rowSums(SetUp[, -1:-2])==0), ] # remove control plots (wild)

SetUp[,-1:-2] <- SetUp[,-1:-2]*1/rowSums(SetUp[,-1:-2]) # convert from present/absent to proportions (All sown evenly)

# group species by functional traits
SetUp$NLF <- rowSums(SetUp[, c(3, 7, 8, 17)]) # non-leguminous forbs
SetUp$C3G <- rowSums(SetUp[, c(4, 10, 11, 15)]) # C3 grass
SetUp$LF <- rowSums(SetUp[, c(5, 12, 13, 14)]) # leguminous forbs
SetUp$C4G <- rowSums(SetUp[, c(6, 9, 16, 18)]) # C4 grass


BioCON <- BioCON[which(BioCON$`Experiment`=='M'),] # select experiment type M

BioCON <- BioCON[which(BioCON$`Sampling..` %% 2 == 0), ] # select one reading per year

# select treatment levels (all ambient/control)
BioCON <- BioCON[which(BioCON$`CO2.Treatment` == "Camb"), ]
BioCON <- BioCON[which(BioCON$`Nitrogen.Treatment` == "Namb"), ]
BioCON <- BioCON[which(BioCON$`Water.Treatment.` == " "), ]
BioCON <- BioCON[which(BioCON$`Temp.Treatment.` == " "), ]

BioCON <- BioCON[which(BioCON$`Sampling..` == 2), ] # select first year of data

# remove unneeded columns
BioCON <- BioCON[, !(names(BioCON) %in% 
                c("Date", "monospecies", "Monogroup", "Experiment", "Sampling..", 
                  "CO2.Treatment", "Nitrogen.Treatment", "Water.Treatment.", "Temp.Treatment."))]


# create columns for achieved biomass of each functional groups
SetUp$NLF_Ach <- 0
SetUp$C3G_Ach <- 0
SetUp$LF_Ach <- 0
SetUp$C4G_Ach <- 0
SetUp$Weeds_Ach <- 0


for(i in 1:nrow(BioCON)) #Remove whitespace from end of response names
{
  if(substr(BioCON$Species[i], nchar(BioCON$Species[i]), nchar(BioCON$Species[i])) == ' ')
  {
    BioCON$Species[i] <- substr(BioCON$Species[i], 1, nchar(BioCON$Species[i])-1)
  }
}

# add biomass of individual species to their functional group columns 
for(i in 1:nrow(BioCON))
{
  plotNb <- BioCON$Plot[i] 
  currentSp <- BioCON$Species[i]
  
  if(currentSp %in% c("Achillea millefolium", "Anemone cylindrica", "Asclepias tuberosa", "Solidago rigida"))
  {
    SetUp[which(SetUp$Plot == plotNb), "NLF_Ach"] <- SetUp[which(SetUp$Plot == plotNb), "NLF_Ach"] + BioCON$Aboveground.Biomass..g.m.2.[i]
  }
  else if(currentSp %in% c("Agropyron repens", "Bromus inermis", "Koeleria cristata", "Poa pratensis"))
  {
    SetUp[which(SetUp$Plot == plotNb), "C3G_Ach"] <- SetUp[which(SetUp$Plot == plotNb), "C3G_Ach"] + BioCON$Aboveground.Biomass..g.m.2.[i]
  }
  else if(currentSp %in% c("Amorpha canescens", "Lespedeza capitata", "Lupinus perennis", "Petalostemum villosum"))
  {
    SetUp[which(SetUp$Plot == plotNb), "LF_Ach"] <- SetUp[which(SetUp$Plot == plotNb), "LF_Ach"] + BioCON$Aboveground.Biomass..g.m.2.[i]
  }
  else if(currentSp %in% c("Andropogon gerardi", "Bouteloua gracilis", "Schizachyrium scoparium", "Sorghastrum nutans"))
  {
    SetUp[which(SetUp$Plot == plotNb), "C4G_Ach"] <- SetUp[which(SetUp$Plot == plotNb), "C4G_Ach"] + BioCON$Aboveground.Biomass..g.m.2.[i]
  }
  else
  {
    SetUp[which(SetUp$Plot == plotNb), "Weeds_Ach"] <- SetUp[which(SetUp$Plot == plotNb), "Weeds_Ach"] + BioCON$Aboveground.Biomass..g.m.2.[i]
  }
}

# merge datasets
BioCON <- BioCON[!duplicated(BioCON[-5:-6]), -5:-6]
BioCON <- merge(BioCON, SetUp[, -2:-18], by = "Plot")

# melt data into stacked format (one row per response reading)
predictors <- setdiff(colnames(BioCON), c("NLF_Ach", "C3G_Ach", "LF_Ach", "C4G_Ach", "Weeds_Ach"))
BioCON <- reshape2::melt(BioCON, id.vars = predictors, variable.name = "YieldSpecies", value.name = "YieldValue")


# use the DImodels package to add DI model interaction columns
BioCON <- cbind(BioCON, DImodels::DI_data(prop = c("NLF", "C3G", "LF", "C4G"), what = c("AV", "ADD"), data = BioCON))

# create a presence/absence column for the functional group response
for(i in 1:nrow(BioCON))
{
  if((BioCON$YieldSpecies[i] == "NLF_Ach") & (BioCON$NLF[i] == 0))
  {
    BioCON$Present[i] <- 0
  }
  else if((BioCON$YieldSpecies[i] == "C3G_Ach") & (BioCON$C3G[i] == 0))
  {
    BioCON$Present[i] <- 0
  }
  else if((BioCON$YieldSpecies[i] == "LF_Ach") & (BioCON$LF[i] == 0))
  {
    BioCON$Present[i] <- 0
  }
  else if((BioCON$YieldSpecies[i] == "C4G_Ach") & (BioCON$C4G[i] == 0))
  {
    BioCON$Present[i] <- 0
  }
  else #Mixture or Weed Response
  {
    BioCON$Present[i] <- 1
  }
}

################################################################################################################################################################

##Continuous Achieved DI models ################################################################################################################################

IDmodel <- gls(YieldValue ~ 0 + Present:(YieldSpecies:(NLF + C3G + LF + C4G)), 
               weights = nlme::varIdent(form = ~ 0 | YieldSpecies), 
               correlation = nlme::corSymm(form = as.formula(~ 0 | Plot)), 
               method = "ML", control=nlme::glsControl(msMaxIter = 50000, opt = "optim"), 
               data=BioCON)

AVmodel <- gls(YieldValue ~ 0 + Present:(YieldSpecies:(NLF + C3G + LF + C4G + AV)), 
               weights = nlme::varIdent(form = ~ 0 | YieldSpecies), 
               correlation = nlme::corSymm(form = as.formula(~ 0 | Plot)), 
               method = "ML", control=nlme::glsControl(msMaxIter = 50000, opt = "optim"), 
               data=BioCON)

ADDmodel <- gls(YieldValue ~ 0 + Present:(YieldSpecies:(NLF + C3G + LF + C4G + ADD.NLF_add + ADD.C3G_add + ADD.LF_add + ADD.C4G_add)), 
                weights = nlme::varIdent(form = ~ 0 | YieldSpecies), 
                correlation = nlme::corSymm(form = as.formula(~ 0 | Plot)), 
                method = "ML", control=nlme::glsControl(msMaxIter = 50000, opt = "optim"), 
                data=BioCON)

FULLmodel <- gls(YieldValue ~ 0 + Present:(YieldSpecies:(NLF + C3G + LF + C4G + (NLF + C3G + LF + C4G)^2)), 
                 weights = nlme::varIdent(form = ~ 0 | YieldSpecies), 
                 correlation = nlme::corSymm(form = as.formula(~ 0 | Plot)), 
                 method = "ML", control=nlme::glsControl(msMaxIter = 50000, opt = "optim"), 
                 data=BioCON)

################################################################################################################################################################

##Model Selection ##############################################################################################################################################

# LRT
anova(IDmodel, AVmodel)
anova(AVmodel, ADDmodel)
anova(AVmodel, FULLmodel) # AVmodel fits the data best

# AIC
AIC(IDmodel); AIC(AVmodel); AIC(ADDmodel); AIC(FULLmodel) # AVmodel fits the data best

# BIC
BIC(IDmodel); BIC(AVmodel); BIC(ADDmodel); BIC(FULLmodel) # AVmodel fits the data best

################################################################################################################################################################

##Final Model ##################################################################################################################################################


# refit the final model using REML to better estimate variances
AVmodel <- gls(YieldValue ~ 0 + Present:(YieldSpecies:(NLF + C3G + LF + C4G + AV)), 
              weights = nlme::varIdent(form = ~ 0 | YieldSpecies), 
              correlation = nlme::corSymm(form = as.formula(~ 0 | Plot)), 
              method = "REML", control=nlme::glsControl(msMaxIter = 50000, opt = "optim"), 
              data=BioCON)

summary(AVmodel)

################################################################################################################################################################

##Plotting #####################################################################################################################################################

# select communities to predict from
plots <- c(67, 70, 77, 100, 109, 121, 309)
communities <- c("0.25:0.25:0.25:0.25", "0:0:0:1", "1:0:0:0", "0:0:1:0", "0:1:0:0", "0:0.5:0:0.5", "0.5:0.25:0.25:0")
richness <- c(4, 1, 1, 1, 1, 2, 3)

predictData <- BioCON[which(BioCON$Plot %in% plots), -12:-15]

# predict from selected communities
predictions <- predict(object = AVmodel, newdata = predictData)
intervals(AVmodel)

# bind information needed for plotting
predictions <- cbind(BioCON[which(BioCON$Plot %in% plots), -12:-15], predictions)
predictions <- cbind(predictions, Community = rep(communities, times = 5))
predictions <- cbind(predictions, richness = rep(richness, times = 5))

predictions$`Community` <- factor(predictions$`Community`, levels = c("1:0:0:0", "0:1:0:0", "0:0:1:0", "0:0:0:1",  "0:0.5:0:0.5", "0.5:0.25:0.25:0", "0.25:0.25:0.25:0.25"))


ggplot2::ggplot(predictions, aes(fill=YieldSpecies, y=predictions, x=Community)) +
         geom_bar(position="stack", stat="identity") +
         scale_fill_manual(values = c("#ffcb66", "#669aff", "#cb66ff", "#9aff66", "#000000"), breaks=c("NLF_Ach", "C3G_Ach", "LF_Ach",  "C4G_Ach", "Weeds_Ach")) +
         facet_grid(~richness, scales = "free_x", space = "free_x") + # add , switch="both" to move richness group name labels to bottom
         expand_limits(x = 2.1) +
         theme(axis.text = element_text(size=20, angle = 20),
               axis.title = element_text(size=24),
               legend.position="top", #remove to put to default right-hand side
               legend.text = element_text(size=20),
               legend.title = element_text(size=20),
               strip.text.x = element_text(size=20),
               panel.border = element_blank(), #white background with no panels/lines
               panel.grid.major = element_blank(), #^
               panel.grid.minor = element_blank(), #^
               panel.background = element_rect(fill = "white"), #^
               panel.spacing = unit(0.25, "cm")) +
         labs(y = "Predicted Yield")

################################################################################################################################################################
################################################################################################################################################################




##ZADR ACHIEVED DI MODEL
##Data Preparation #############################################################################################################################################

setwd("C:/Users/Administrator/OneDrive - TCDUD.onmicrosoft.com/Documents/PhD Stuff/Tasks/Task 6 - RPackage/COST database")

AgroDiv <- read.csv("biomass.csv") # biomass response data

AgroDiv <- AgroDiv[which(AgroDiv$SITE == 1) ,] # Belgium site
AgroDiv <- AgroDiv[which(AgroDiv$YEARN == 1) ,] # Year 2003 (first harvest year)

AgroDiv <- AgroDiv %>% group_by(PLOT, G1, G2, L1, L2, TREAT, DENS,) %>%
                       summarise(G1_Y = sum(G1_Y),
                                 G2_Y = sum(G2_Y),
                                 L1_Y = sum(L1_Y),
                                 L2_Y = sum(L2_Y),
                                 WEED_Y = sum(WEED_Y),
                                 HARV_YIELD = sum(HARV_YIELD)
                                 ) %>% as.data.frame()

AgroDiv <- cbind(AgroDiv, DI_data(prop = 2:5, FG = c("G", "G", "L", "L"), what = c("AV", "ADD", "FG", "FULL"), data = AgroDiv))

AgroDiv$G1_Yp <- AgroDiv$G1_Y / AgroDiv$HARV_YIELD
AgroDiv$G2_Yp <- AgroDiv$G2_Y / AgroDiv$HARV_YIELD
AgroDiv$L1_Yp <- AgroDiv$L1_Y / AgroDiv$HARV_YIELD
AgroDiv$L2_Yp <- AgroDiv$L2_Y / AgroDiv$HARV_YIELD
AgroDiv$WEED_Yp <- AgroDiv$WEED_Y / AgroDiv$HARV_YIELD


################################################################################################################################################################

##ZADR Achieved DI models #####################################################################################################################################

Yp <- as.matrix(AgroDiv[, 28:32])
X_ID <- as.matrix(AgroDiv[, 2:5])

IDmodel <- Compositional::zadr2(y = Yp, x = X_ID, con = FALSE)


Yp <- as.matrix(AgroDiv[, 28:32])
X_AV <- as.matrix(AgroDiv[, c(2:5, 14)])

AVmodel <- Compositional::zadr2(y = Yp, x = X_AV, con = FALSE)


X_FG <- as.matrix(AgroDiv[, c(2:5, 15:17)])

FGmodel <- Compositional::zadr2(y = Yp, x = X_FG, con = FALSE)


X_ADD <- as.matrix(AgroDiv[, c(2:5, 18:21)])

ADDmodel <- Compositional::zadr2(y = Yp, x = X_ADD, con = FALSE)


X_FULL <- as.matrix(AgroDiv[, c(2:5, 22:27)])

FULLmodel <- Compositional::zadr2(y = Yp, x = X_FULL, con = FALSE)



Tmodel <- DImodels::DI(prop = 2:5, y = 13, DImodel = "AV", data = AgroDiv[, -14:-32])


################################################################################################################################################################

##Model Selection #############################################################################################################################################

#ID vs AV (df = 1)
pchisq((-2 * (IDmodel$loglik - AVmodel$loglik)), df = 1, lower.tail = FALSE)
#AV model Better

#AV vs ADD (df = 3)
pchisq((-2 * (AVmodel$loglik - ADDmodel$loglik)), df = 3, lower.tail = FALSE)
#AV model Better

#AV vs FG (df = 2)
pchisq((-2 * (AVmodel$loglik - FGmodel$loglik)), df = 2, lower.tail = FALSE)
#AV model Better

#AV vs FULL (df = 5)
pchisq((-2 * (AVmodel$loglik - FULLmodel$loglik)), df = 5, lower.tail = FALSE)
#AV model Better


#AIC
2 * (length(IDmodel$be) + 1) - 2 * IDmodel$loglik

2 * (length(AVmodel$be) + 1) - 2 * AVmodel$loglik

2 * (length(ADDmodel$be) + 1) - 2 * ADDmodel$loglik

2 * (length(FGmodel$be) + 1) - 2 * FGmodel$loglik

2 * (length(FULLmodel$be) + 1) - 2 * FULLmodel$loglik


#BIC
log(nrow(X_ID)) * (length(IDmodel$be) + 1) - 2 * IDmodel$loglik

log(nrow(X_AV)) * (length(AVmodel$be) + 1) - 2 * AVmodel$loglik

log(nrow(X_ADD)) * (length(ADDmodel$be) + 1) - 2 * ADDmodel$loglik

log(nrow(X_FG)) * (length(FGmodel$be) + 1) - 2 * FGmodel$loglik

log(nrow(X_FULL)) * (length(FULLmodel$be) + 1) - 2 * FULLmodel$loglik


################################################################################################################################################################

##Results ######################################################################################################################################################

AVmodel$be

AVmodel$sigma

summary(Tmodel)

################################################################################################################################################################

##Plots ########################################################################################################################################################


################################################################################################################################################################

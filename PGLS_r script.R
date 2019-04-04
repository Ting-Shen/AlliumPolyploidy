####################################################################################################
#################################### PGLS Tests By Ting-Shen Han ###################################

#pgls for continuous traits

#Remove all data

rm(list=ls())

#Set your working directory where all data files are located

setwd("...")

#Load relevant libraries or packages

library(ape)
library(caper)
library(geiger)
library(nlme)
library(phytools)

#load character matrix and tree

whetPoly<-read.csv("pglsdata.csv", row.names=1,header = T)
tree1<-read.tree("mytree.tre")

# drop.tip of species in tree and tips with no data

whetPoly
whetPoly.v<-whetPoly[,1]
names(whetPoly.v)<-row.names(whetPoly)
x<-name.check(tree1, whetPoly)
tree2<-drop.tip(tree1,c(x$tree_not_data))

#Comparative dataset creation

pglsdata <- comparative.data (tree2, whetPoly, names.col="Species.1", vcv=TRUE, vcv.dim = 3)
pglsdata

#model selection for the best transformation structure of the covariance matrix
#response variable ~ explanatory variable

mod.l <- pgls(ClimatePC2 ~ PolyploidPC2, pglsdata, lambda = 'ML', param.CI = 0.95)
mod.d <- pgls(ClimatePC2 ~ PolyploidPC2, pglsdata, delta = 'ML', param.CI = 0.95)
mod.k <- pgls(ClimatePC2 ~ PolyploidPC2, pglsdata, kappa = 'ML', param.CI = 0.95)

summary (mod.l)
summary (mod.d)
summary (mod.k)

#pgls

pglsModel <- pgls(ClimatePC2 ~ PolyploidPC2, pglsdata, lambda = 0.154, delta = 3, kappa=0.029, param.CI = 0.95)

summary(pglsModel)

#print the best model

mod1.sum <- summary(mod1)
print(mod1.sum)

shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
mod1 <- pgls(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird, lambda='ML')
mod2 <- pgls(log(Egg.Mass) ~ log(M.Mass), data=shorebird, lambda='ML', delta='ML')
mod3 <- pgls(log(Egg.Mass) ~ log(M.Mass), data=shorebird, lambda = 1.0, kappa = 1.0,  delta= 1.0, param.CI = 0.95, control = list(fnscale=-1), bounds = NULL)

#run phylogenetic generalized least squares in nlme

pglsModel <- gls(PolyploidPC2 ~ SoilPC2, correlation = corBrownian(phy = tree3), data = whetPoly, method = "ML")
summary(pglsModel)

#Remove all data

rm(list=ls())

#Set your working directory where all data files are located

setwd("...")

#Load relevant libraries or install packages:

library(ape)
library(geiger)
library(mvtnorm)

#Load trait data and phylogeny. In this case 'fruit_length' refers to our continuous fruit trait data.

#trait

trait0<-read.csv("mydata.CSV", row.names=1, header = FALSE)
trait.v<-trait0[,1]
names(trait.v)<-row.names(trait0)

#phylogeny in nexus format

tree2<-read.tree("mytree.tre")

#Remove tips from tree if absent from trait data file

x<-name.check(tree2, trait0)
tree3<-drop.tip(tree2,c(x$tree_not_data))

#Data names

trait <- trait.v
phy <-tree3


#simply download essim.R from the ./R/ directory of this repository, then load it in R:
#Run it using the command

#es_pgls

source(".../espgls.R")
espgls(phy, trait)

#nd_pgls

source(".../ndpgls.R")
ndpgls(phy, trait)

#tb_pgls

source(".../tbpgls.R")
tbpgls(phy, trait)

#nd_sim

source(".../ndsim.R")
ndsim(phy, trait, nsim = 1000)

#tb_sim

source(".../tbsim.R")
tbsim(phy, trait, nsim = 1000)

#essim_pearson

source("essim.R")
essim(phy, trait, nsim = 1000)

#essim_spearman

source("essim_spearman.R")
essim_spearman(phy, trait, nsim = 1000)

#essim_kendal

source(".../essim_kendall.R")
essim_kendall(phy, trait, nsim = 1000)

#Run it using the command

essim(phy, trait, nsim = 1000)
espgls(phy, trait)
essim_kendall(phy, trait, nsim = 1000)

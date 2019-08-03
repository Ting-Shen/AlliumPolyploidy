####################################################################################################
################################## MuHiSSE Tests By Ting-Shen Han ##################################

rm(list=ls())

# set the working directory

setwd("...")

# load packages

library(hisse)
library(diversitree)
library(geiger)
library(phytools)

# load phylogeny and trait data

mytree <- read.tree("mytree.tre")
mytrait <- read.csv("mytrait.csv")

# make a new dateframe of the states and then convert the first column to characters rather than factors (Harrington & Reeder, 2017)

statesALL<-mytrait
statesALL[,1]<-as.character(mytrait[,1])
rownames(statesALL)<-statesALL[,1]

# remove data for taxa not in the tree and tips with no data (Harrington & Reeder, 2017)

check<-name.check(mytree, statesALL)
states<-(statesALL[!rownames(statesALL) %in% check$data_not_tree,])
states

phy<-drop.tip(mytree, tip=check$tree_not_data)
name.check(phy, states)

# sampling fraction for each state
# calculated from imputation result of phylopars

sampfrac <- c(xxx, xxx)

#### set up and execute MuHiSSE models, and estimate the likeliest states

# MuHiSSE model with no dual and qPD=0, a MuSSE-like full model (MuL.qPD0)

trans.rate.nodual <- TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.nodual

trans.rate.irrev1 <- ParDrop(trans.rate.nodual, c(2,4,5,7))
trans.rate.irrev1

MuL.qPD0 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(1,2,3,4), eps=c(1,2,3,4), hidden.states=FALSE, trans.rate=trans.rate.irrev1)


# MuHiSSE model without dual transitions between multistate traits, this is a null model, MuSSE-like.nodual.null (Mul.nodual.null)

trans.rate.nodual <- TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.nodual


# MuHiSSE model with a true MuSSE model, this is a MuSSE-like.nodual.full model 01 (Mul.nodual.full01)

Mul.nodual.full01 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(1,2,3,4), eps=c(1,1,1,1), hidden.states=FALSE, trans.rate=trans.rate.nodual)


# MuHiSSE model with a true MuSSE model, this is a MuSSE-like.nodual.full model 02 (Mul.nodual.full02)

Mul.nodual.full02 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(1,2,3,4), eps=c(1,2,3,4), hidden.states=FALSE, trans.rate=trans.rate.nodual)


# MuHiSSE model with character-dependent full model 01, with 2 hidden states (MuH.nodual.full01)

trans.rate.nodual.h1 <- TransMatMakerMuHiSSE(hidden.traits=1)
trans.rate.nodual.h1

MuH.nodual.full01 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(1,2,3,4,5,6,7,8), eps=c(1,1,1,1,1,1,1,1), hidden.states=TRUE, trans.rate=trans.rate.nodual.h1)


# MuHiSSE model with character-dependent full model 02, with 2 hidden states (MuH.nodual.full02)

MuH.nodual.full02 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(1,2,3,4,5,6,7,8), eps=c(1,2,3,4,5,6,7,8), hidden.states=TRUE, trans.rate=trans.rate.nodual.h1)


## Allow for dual transitions, e.g., 00->11 and 11->00 transitions, so repeat the above models once more

trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0, include.diagonals = TRUE)
trans.rate


# MuHiSSE model with a true MuSSE model, this is a MuSSE-like.full model 01 (Mul.full01)

Mul.full01 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(1,2,3,4), eps=c(1,1,1,1), hidden.states=FALSE, trans.rate=trans.rate)


# MuHiSSE model with a true MuSSE model, this is a MuSSE-like.full model 02 (Mul.full02)

Mul.full02 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(1,2,3,4), eps=c(1,2,3,4), hidden.states=FALSE, trans.rate=trans.rate)


# MuHiSSE model with character-dependent full model 01, with 2 hidden states (MuH.full01)

trans.rate.h1 <- TransMatMakerMuHiSSE(hidden.traits=1, include.diagonals = TRUE)
trans.rate.h1

MuH.full01 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(1,2,3,4,5,6,7,8), eps=c(1,1,1,1,1,1,1,1), hidden.states=TRUE, trans.rate=trans.rate.h1)


# MuHiSSE model with character-dependent full model 02, with 2 hidden states (MuH.nodual.full02)

MuH.full02 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(1,2,3,4,5,6,7,8), eps=c(1,2,3,4,5,6,7,8), hidden.states=TRUE, trans.rate=trans.rate.h1)


# MuHiSSE model with qPD=0, full model 02 (MuH.qPD0)

trans.rate.irrev2 <- ParDrop(trans.rate.h1, c(2,3,5,6,7,8,10,11,14,15,17,18,19,20,22,23))
trans.rate.irrev2

MuH.qPD0 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(1,2,3,4,5,6,7,8), eps=c(1,2,3,4,5,6,7,8), hidden.states=TRUE, trans.rate=trans.rate.irrev2)


# MuHiSSE model with character-independent model, with 2 hidden states (MuH.CID2)

MuH.CID2 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(rep(1,4),rep(2,4)), eps=c(rep(1,4),rep(2,4)), hidden.states=TRUE, trans.rate=trans.rate.h1)


# MuHiSSE model with character-independent model, with 3 hidden states (MuH.CID3)

trans.rate.h2 <- TransMatMakerMuHiSSE(hidden.traits=2, include.diagonals = TRUE)
trans.rate.h2

MuH.CID3 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(rep(1,4),rep(2,4),rep(3,4)), eps=c(rep(1,4),rep(2,4),rep(3,4)), hidden.states=TRUE, trans.rate=trans.rate.h2)


# MuHiSSE model with character-independent model, with 4 hidden states (MuH.CID4)

trans.rate.h3 <- TransMatMakerMuHiSSE(hidden.traits=3, include.diagonals = TRUE)
trans.rate.h3

MuH.CID4 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4)), eps=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4)), hidden.states=TRUE, trans.rate=trans.rate.h3)


# MuHiSSE model with character-independent model, with 5 hidden states (MuH.CID5)

trans.rate.h4 <- TransMatMakerMuHiSSE(hidden.traits=4, include.diagonals = TRUE)
trans.rate.h4

MuH.CID5 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4)), eps=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4)), hidden.states=TRUE, trans.rate=trans.rate.h4)


# MuHiSSE model with character-independent model, with 6 hidden states (MuH.CID6)

trans.rate.h5 <- TransMatMakerMuHiSSE(hidden.traits=5, include.diagonals = TRUE)
trans.rate.h5

MuH.CID6 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4)), eps=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4)), hidden.states=TRUE, trans.rate=trans.rate.h5)


# MuHiSSE model with character-independent model, with 7 hidden states (MuH.CID7)

trans.rate.h6 <- TransMatMakerMuHiSSE(hidden.traits=6, include.diagonals = TRUE)
trans.rate.h6

MuH.CID7 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4)), eps=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4)), hidden.states=TRUE, trans.rate=trans.rate.h6)


# MuHiSSE model with character-independent model, with 8 hidden states (MuH.CID8)

trans.rate.h7 <- TransMatMakerMuHiSSE(hidden.traits=7, include.diagonals = TRUE)
trans.rate.h7

MuH.CID8 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4),rep(8,4)), eps=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4),rep(8,4)), hidden.states=TRUE, trans.rate=trans.rate.h7)


# MuHiSSE model with character-independent model, with no dual transition, and 2 hidden states (MuH.nodual.CID2)

trans.rate.nodual.h1 <- TransMatMakerMuHiSSE(hidden.traits=1)
trans.rate.nodual.h1

MuH.nodual.CID2 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(rep(1,4),rep(2,4)), eps=c(rep(1,4),rep(2,4)), hidden.states=TRUE, trans.rate=trans.rate.nodual.h1)


# MuHiSSE model with character-independent model, with no dual transition, and 4 hidden states (MuH.nodual.CID4)

trans.rate.nodual.h3 <- TransMatMakerMuHiSSE(hidden.traits=3)
trans.rate.nodual.h3

MuH.nodual.CID4 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4)), eps=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4)), hidden.states=TRUE, trans.rate=trans.rate.nodual.h3)


# MuHiSSE model with character-independent model, with no dual transition, and 8 hidden states (MuH.nodual.CID8)

trans.rate.nodual.h7 <- TransMatMakerMuHiSSE(hidden.traits=7)
trans.rate.nodual.h7

MuH.nodual.CID8 <- MuHiSSE(phy=phy, data=states, f=sampfrac, turnover=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4),rep(8,4)), eps=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4),rep(8,4)), hidden.states=TRUE, trans.rate=trans.rate.nodual.h7)


# List of all recons which can be plotted

recon_list<-list(MuL.qPD0, Mul.nodual.full01, Mul.nodual.full02, MuH.nodual.full01, MuH.nodual.full02, Mul.full01, Mul.full02, MuH.full01, MuH.full02, MuH.qPD0, MuH.CID2, MuH.CID3, MuH.CID4, MuH.CID5, MuH.CID6, MuH.CID7, MuH.CID8, MuH.nodual.CID2, MuH.nodual.CID4, MuH.nodual.CID8)
recon_list

# Adaptive Sampling of the Likelihood Surface under the best fitted MuHiSSE (e.g., MuH.full01)

MuH.full01.CI <- SupportRegionMuHiSSE(muhisse.obj = MuH.full01, n.points = 10000)
write.csv(MuH.full01$all.points, "MuHiSSE_MuH.full01.CI.csv")

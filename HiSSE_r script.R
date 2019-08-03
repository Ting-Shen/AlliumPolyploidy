####################################################################################################
################################### HiSSE Tests By Ting-Shen Han ###################################

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

#### set up and execute HiSSE models, and estimate the likeliest states
## BiSSE-like models (no hidden states)

# BiSSE-like unequal turnover rate and extinction rate (BiSSE-like full model, BLFull)

BLFull<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=TransMatMaker())


# BiSSE-like equal turnover rate but unequal extinction rate (BiSSE-like equal.t model, BLEqual.t)

BLEqual.t<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,2,0,0), trans.rate=TransMatMaker())


# BiSSE-like unequal turnover rate but equal extinction rate (BiSSE-like equal.e model, BLEqual.e)

BLEqual.e<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=TransMatMaker())


# BiSSE-like equal transition (BiSSE-like equal.q model, BLEqual.q)

trans.rates<-TransMatMaker()
trans.rates
trans.rates.allequal<-ParEqual(trans.rates, c(1,2))
trans.rates.allequal

BLEqual.q<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.allequal)


# BiSSE-like non transition from polyploidy to diploidy (BiSSE-like qPD0 model, BL.qPD0)

trans.rates<-TransMatMaker()
trans.rates
trans.rates.irrev<-ParDrop(trans.rates, 1)
trans.rates.irrev

BL.qPD0 <- hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.irrev)


# BiSSE-like equal turnover rate and extinction rate (BiSSE-like equal.te model, BLEqual.te)

BLEqual.te<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.allequal)


# BiSSE-like equal turnover rate and transition rate (BiSSE-like equal.tq model, BLEqual.tq)

trans.rates<-TransMatMaker()
trans.rates
trans.rates.allequal<-ParEqual(trans.rates, c(1,2))
trans.rates.allequal

BLEqual.tq<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.allequal)


# BiSSE-like equal extinction rate and transition rate (BiSSE-like equal.eq model, BLEqual.eq)

trans.rates<-TransMatMaker()
trans.rates
trans.rates.allequal<-ParEqual(trans.rates, c(1,2))
trans.rates.allequal

BLEqual.eq<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.allequal)


# BiSSE-like equal turnover rate and non transition from polyploidy to diploidy (BiSSE-like equal.t, qPD0 model, BLEqual.t.qPD0)

trans.rates<-TransMatMaker()
trans.rates
trans.rates.irrev<-ParDrop(trans.rates, 1)
trans.rates.irrev

BLEqual.t.qPD0<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.irrev)


# BiSSE-like equal extinction rate and non transition from polyploidy to diploidy (BiSSE-like equal.e, qPD0 model, BLEqual.e.qPD0)

trans.rates<-TransMatMaker()
trans.rates
trans.rates.irrev<-ParDrop(trans.rates, 1)
trans.rates.irrev

BLEqual.e.qPD0<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.irrev)


# BiSSE-like equal turnover rate, extinction rate and transition rate (BiSSE-like equal.teq model, BLEqual.teq)

BLEqual.teq<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0), trans.rate=TransMatMaker())


# BiSSE-like equal turnover rate and equal extinction rate, but non transition from polyploidy to diploidy (BiSSE-like equal.te, qPD0 model, BLEqual.te.qPD0)

trans.rates<-TransMatMaker()
trans.rates
trans.rates.irrev<-ParDrop(trans.rates, 1)
trans.rates.irrev

BLEqual.te.qPD0<-hisse(phy, states, f=sampfrac, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.irrev)


## HiSSE models (with hidden states)
# HiSSE full model (HFull)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates

HFull<-hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates)


# HiSSE model without dual transitions between both the observed trait and the hidden trait (HiSSE nodual, H.nodual)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates
trans.rates.nodual <- ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual

H.nodual<-hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual)


# HiSSE model with all equal transitions (HiSSE equal.q, H.Equal.q)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates
trans.rates.allequal <- ParEqual(trans.rates, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,11,1,12))
trans.rates.allequal

H.Equal.q <- hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.allequal)


# HiSSE model with equal transitions for nodual ones (HiSSE equal.q/nodual, H.nodual.Equal.q)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates
trans.rates.nodual.allequal <- ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

H.nodual.Equal.q <- hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual.allequal)


# HiSSE model with non transition from polyploidy to diploidy (HiSSE qPD0, H.qPD0)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates
trans.rates.irrev2 <- ParDrop(trans.rates, c(1,3,8,9))
trans.rates.irrev2

H.qPD0 <- hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.irrev2)


# HiSSE model with non dual transition and non transition from polyploidy to diploidy (HiSSE qPD0/nodual, H.nodual.qPD0)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates
trans.rates.nodual <- ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
trans.rates.irrev3 <- ParDrop(trans.rates.nodual, c(1,0,6))
trans.rates.irrev3

H.nodual.qPD0 <- hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.irrev3)


# HiSSE model with all equal transitions and non transition from polyploidy to diploidy (HiSSE euqal.q/qPD0, H.Equal.q.qPD0)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates
trans.rates.irrev4 <- ParDrop(trans.rates, c(1,3,8,9))
trans.rates.irrev4
trans.rates.allequal2 <- ParEqual(trans.rates.irrev4, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.allequal2

H.Equal.q.qPD0 <- hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.allequal2)


# HiSSE model with equal transitions for nodual ones and non transition from polyploidy to diploidy (HiSSE equal.q/qPD0/nodual, H.nodual.Equal.q.qPD0)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates
trans.rates.nodual2 <- ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual2
trans.rates.irrev5 <- ParDrop(trans.rates.nodual2, c(1,0,6))
trans.rates.irrev5
trans.rates.nodual.allequal2 <- ParEqual(trans.rates.irrev5, c(1,2,1,3,1,4,1,5))
trans.rates.nodual.allequal2

H.nodual.Equal.q.qPD0 <- hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual.allequal2)


# HiSSE model with equal turnover rate (HiSSE equal.t, H.Equal.t)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates

H.Equal.t<-hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,1,1,1), eps.anc=c(1,2,3,4), trans.rate=trans.rates)


# HiSSE model with equal extinction rate (HiSSE equal.e, H.Equal.e)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates

H.Equal.e<-hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates)


# CID2
# HiSSE null-two model (HiSSE null-2,HNull2)

trans.rates.nodual[!is.na(trans.rates.nodual) & !trans.rates.nodual == 0] = 1

HNull2 <- hisse(phy, states, f=sampfrac, hidden.states = TRUE, turnover.anc=c(1,1,2,2),eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual)


# HiSSE null-two model with non transition from polyploidy to diploidy (HiSSE null-2/qPD0, HNull2.qPD0)

HNull2.qPD0 <- hisse(phy, states, f=sampfrac, hidden.states = TRUE, turnover.anc=c(1,1,2,2),eps.anc=c(1,1,2,2), trans.rate=trans.rates.irrev2)


# HiSSE null-two model with non dual transition (HiSSE null-2/nodual, HNull2.nondual)

HNull2.nondual <- hisse(phy, states, f=sampfrac, hidden.states = TRUE, turnover.anc=c(1,1,2,2),eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual)


# CID3
# HiSSE null-three model, HNull3

HNull3 <- hisse.null4(phy, states, f=sampfrac, trans.type = "three.rate")


# CID4 
# HiSSE null-four model, HNull4

HNull4 <- hisse.null4(phy, states, f=sampfrac, trans.type = "equal")


# HiSSE model with 3 rates (HiSSE 3 rates, H.3rates) 

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates
trans.rates.nodual <- ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
trans.rates.nodual.threerates <- trans.rates.nodual

to.change <- cbind(c(1,3), c(2,4))
trans.rates.nodual.threerates[to.change] = 1

to.change <- cbind(c(2,4), c(1,3))
trans.rates.nodual.threerates[to.change] = 2

to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.nodual.threerates[to.change] = 3
trans.rates.nodual.threerates

H.3rates<-hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.nodual.threerates)


# HiSSE model with two rate irreversible (HiSSE 2 rate/irreversible, H.2rates.irrev)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates
trans.rates.nodual <- ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
trans.rates.two.rates<- trans.rates.nodual
to.change <- cbind(c(1,3), c(2,4))
trans.rates.two.rates[to.change] = 1

to.change <- cbind(c(2,4), c(1,3))
trans.rates.two.rates[to.change] = 0

to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.two.rates[to.change] = 2
trans.rates.two.rates

H.2rates.irrev <-hisse(phy, states,  f=sampfrac, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.two.rates)


# HiSSE null-two model with three rates (HiSSE null-2/3 rates, HNull2.3rates)

HNull2.3rates<- hisse(phy, states, f=sampfrac, hidden.states = TRUE, turnover.anc=c(1,1,2,2),eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual.threerates)


# List of all recons which can be plotted

recon_list<-list(BiSSE_BD_states, null.4_9.rate_states, null.4_3.rate_states, BiSSE_states,HiSSE_states,HiSSE_states_irrev,BiSSE_states_Irrev,null.4_Eq.rates_states,HiSSE_EqDiv_states,hisse_res_HS_all_equal_States,null_two_3_rate_States,hisse_res_HS_3_rate_States,hisse_res_HS_all_equal_IRREV_States,hisse_res_HS_2_rate_IRREV_States,null_two_8_rates_States,null_two_Irrev_6_rates_States)

plot.hisse.states(H.nodual.qPD0.states, rate.param = "net.div", show.tip.label = FALSE, rate.colors=c("white", "red"), state.colors = c("yellow","orange"))


# Obtain an estimate of the confidence intervals for HiSSE model (e.g., H.nodual) using adaptive sampling points

H.nodual.CI <- SupportRegion(hisse.obj = H.nodual, n.points = 10000, output.type = "raw")
write.csv(H.nodual.CI$all.points, "HiSSE_H.nodual.CI.csv")

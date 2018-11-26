## ----load_packages, include=FALSE----------------------------------------
knitr::opts_knit$set(progress = TRUE, verbose = TRUE)
library(hicrep)
data("HiCR1")
data("HiCR2")

## ------------------------------------------------------------------------
dim(HiCR1)
HiCR1[1:10,1:10]

## ---- eval=TRUE----------------------------------------------------------
Pre_HiC <- prep(HiCR1, HiCR2, 1000000, 1, 5000000)
head(Pre_HiC)

## ---- eval=TRUE----------------------------------------------------------
h_hat <- htrain(HiCR1, HiCR2, 1000000, 5000000, 0:2)

h_hat

## ------------------------------------------------------------------------
#check total number of reads before adjustment
sum(HiCR1[,-c(1:3)])

DS_HiCR1 <- depth.adj(HiCR1, 200000, 1000000, out = 0) 

#check total number of reads after adjustment
sum(DS_HiCR1[,-c(1:3)])


## ---- eval=TRUE----------------------------------------------------------
SCC.out = get.scc(Pre_HiC, 1000000, 5000000)

#SCC score
SCC.out[[3]]

#Standard deviation of SCC
SCC.out[[4]]
sessionInfo()


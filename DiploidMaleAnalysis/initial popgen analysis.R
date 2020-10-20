##### Population genetic analyses of Bombus hypnorum worker microsat data
## L P Crowther, September 2016
#  Script performs all popgen analyses for workers
#  Outputs table with HWE etc. for all markers, M ratio figure and Input file for Bottleneck

rm(list=ls()) # clear R

# Load packages
library(adegenet)
library(diveRsity)
library(pegas) # for HWE test as no longer available in adegenet
library(poppr) # for subsetting genind objects more easily than seppop in adegenet
library(hierfstat)
library(PopGenReport)
library(devtools)
library(ggplot2)
library(stringi)
library(dplyr)

#genind.2015 <- read.genetix("data out/workers_2015-GENETIX.gtx")
genind.workers.all<- read.genetix("data out/g.workers_all-GENETIX.gtx")
markers <- read.csv(file = "data out/markerlist.csv")
hwe.all <-hw.test(genind.workers.all)
stat.tab.all <- summary(genind.workers.all)
### many loci not in hwe


## some are decidely worse than others, though. 
## However the homo excess varies.
plot(stat.tab.all$Hobs, col = "black", type = "line")
lines(stat.tab.all$Hexp, col = "red")

# These tests treat both year's samples as a single population
# This is not correct as allele frequencies may have changed due to
# migrants from a population spatially structured at larger scales

genind.2014 <- read.genetix("data out/g.workers.2014-GENETIX.gtx")
genind.2015 <- read.genetix("data out/g.workers_2015-GENETIX.gtx")

hwe.14 <- hw.test(genind.2014)
hwe.15 <- hw.test(genind.2015)
stat.tab.14 <- summary(genind.2014)
stat.tab.15 <- summary(genind.2015)

#plot in 2x2
par(mfrow= c(2,2))

plot(stat.tab.all$Hobs, col = "black", type = "line", main = "all years")
lines(stat.tab.all$Hexp, col = "red")

plot(stat.tab.14$Hobs, col = "black", type = "line", main = "2014")
lines(stat.tab.14$Hexp, col = "red")

plot(stat.tab.15$Hobs, col = "black", type = "line", main = "2015")
lines(stat.tab.15$Hexp, col = "red")

res.tab <- data.frame(markers, 
                      HW_ALL_H0 = logical(length = 20), 
                      p_all = numeric(length =20), 
                      HW_2014_H0 = logical(length = 20),
                      p_14 = numeric(length =20),
                      HW_2015_H0 = logical(length = 20),
                      p_15 = numeric(length =20))

names(res.tab)[1] <- "Marker"

cervus.sum.tab <- read.table("data out/workers.all.cervus.summary.txt", header = TRUE)
# based on HObs the loci are in the same order
# take cols from cervus output
res.tab$k <- cervus.sum.tab$k
res.tab$N <- cervus.sum.tab$N
res.tab$F.NUll <- cervus.sum.tab$F.Null.

res.tab <- res.tab[c(1,9,10,2,3,4,5,6,7,11,8)]
View(res.tab)

# Drop markers where Bonf. corrected HWE Null 
# hyp is rejected in both years OR F(NULL) > 0.1
for(i in 1:20){
  res.tab$p_all[i] <- hwe.all[i, 4]
  res.tab$p_14[i] <- hwe.14[i, 4]
  res.tab$p_15[i] <- hwe.15[i, 4]
  if(!hwe.all[i, 4] > (0.05/60)){
    res.tab$HW_ALL_H0[i] <- FALSE
  }else{
    res.tab$HW_ALL_H0[i] <- TRUE
  }
  if(!hwe.14[i, 4] > (0.05/60)){
    res.tab$HW_2015_H0[i] <- FALSE
  }else{
    res.tab$HW_2014_H0[i] <- TRUE
  }
  if(!hwe.15[i, 4] > (0.05/60)){
    res.tab$HW_2015_H0[i] <- FALSE
  }else{
    res.tab$HW_2015_H0[i] <- TRUE
  }
  if(res.tab$HW_2014_H0[i] == FALSE & res.tab$HW_2015_H0[i] == FALSE){
    res.tab$Action[i] <- "Drop Marker"
  }else{
    res.tab$Action[i] <- "Retain Marker"
  }
  if(res.tab$F.NUll[i] > 0.1){
    res.tab$Action[i] <- "Drop Marker"
  }
  }

res.tab
# make a logical vector of pipeline results for chopping up data later
include <- grepl("Retain", res.tab$Action)

## look at linkage

alleles <- read.csv("data out/g.worker.mat_all.csv")
alleles <- alleles[,3:42]

loci <- alleles2loci(alleles, phased = TRUE)
linkage <- LDscan(loci)
summary(linkage)
# warning, large matrix!
linkage

# BTERN01 has some alleles correlated with those of other markers, 
# but we're already dropping it

## allele frequencies
all.freq.tab <- tab(genind.workers.all, freq = TRUE)

# standard counting function
n_not_na <- function(x) length(which(!is.na(x)))

# summarise by columns
all.freq.n <- apply(all.freq.tab, 2, n_not_na)
all.freq.sum <-  vector(mode = "numeric", length = ncol(all.freq.tab))
for(i in 1: ncol(all.freq.tab)){
  all.freq.sum[i] <- sum(all.freq.tab[,i], na.rm = TRUE)/all.freq.n[i]
  }

# rejoin with allele names
names(all.freq.sum) <- colnames(all.freq.tab) 
write.csv( all.freq.sum, file = "data out/allele frequencies.csv")

# mratio function author = Roman Lustrik, roman.lustrik@biolitika.si, used with thanks.
# https://github.com/romunov/zvau/blob/master/R/Mratio.R

Mratio <- function(x) {
  sapply(alleles(x), FUN = function(m) {
    k <- length(m)
    num.m <- as.numeric(m)
    rn.m <- range(num.m)
    R <- diff(rn.m)
    k/(R + 1)
  })
}


# Get files from Dreier et al., on NERC EIDC
dreier_files <- grep("worker", list.files(path = "./Dreier et al data", full.names = TRUE), value = TRUE)
dreier_data <- list()
for(i in 1: length(dreier_files)){
  dreier_data[[i]] <- read.csv(file = dreier_files[i])
}

# take selected workers
cut_workers <- function(x){
  cut <- x[x$Reduced.data == "X",]
  cut
}
dreier_data <- lapply(dreier_data, cut_workers)
# drop unneeded columns
drop_cols <- function(x) x <- x[, c(3, 6 : ncol(x))]
dreier_data <- lapply(dreier_data, drop_cols)

# calculate M ratio
dreier_loci <- lapply(dreier_data, function(x) alleles2loci(x, rownames = 1))
dreier_genind <- lapply(dreier_loci, loci2genind)
dreier_mrat <- lapply(dreier_genind, Mratio)

# use species codes from Stephanie's file names
names(dreier_mrat) <- lapply(dreier_files, function(x) stri_split_fixed(x, pattern = "_")[[1]][5])

dreier_mrat_means <- lapply(dreier_mrat, mean)

## need to filter B. hypnorum to 1 per colony
clusters <- read.csv("./data out/Bombus hypnorum workers.BestCluster", 
                     sep = "", stringsAsFactors = FALSE,
                     na.strings = c("", "NA"))

select_workers <- function(x){
nests <- unique(x$ClusterIndex)
for(i in 1:length(nests)){
  cworkers <- x[x$ClusterIndex==nests[i],]
  selectworker <- cworkers[sample(1:nrow(cworkers), size = 1),]
  if(i == 1){
    output <- selectworker
  }else{
    output <- rbind(output, selectworker)
  }
}
output
}
workers <- select_workers(clusters)
workers <- workers$OffspringID

# select my workers
g.bhy.selected <- genind.workers.all[which(rownames(genind.workers.all@tab) %in% workers),]
# select my loci
g.bhy.selected <- g.bhy.selected[,loc = include]
  
mls <- Mratio(g.bhy.selected)
#plot them together
mrats <- dreier_mrat
mrats[[6]] <-  mls
names(mrats)[[6]] <- "bhyp"
mrats<- unlist(mrats)
mrats_levels <- lapply(names(mrats), function(x) stri_split_fixed(x, pattern = ".")[[1]][1])
m_ratio_spp <- data.frame(species = unlist(mrats_levels), m = mrats)
bin_names <- c("B. hortorum", "B. hypnorum", "B. lapidarius", "B. pascuorum", "B. ruderatus", "B. terrestris")
levels(m_ratio_spp$species) <- bin_names

# calculatate sample sizes for labels
m_ratio_spp %>% group_by(species) %>% summarise(n = n_distinct(m), mean = mean(m), med = median(m)) -> spp_sum
dreier_genind[[6]] <- g.bhy.selected
n_bees_list <- lapply(dreier_genind,  function(x) n_bees <- nrow(x@tab))
names(n_bees_list) <- lapply(dreier_files, function(x) stri_split_fixed(x, pattern = "_")[[1]][5])
names(n_bees_list) <- c("B. hortorum", "B. lapidarius", "B. pascuorum", "B. ruderatus", "B. terrestris", "B. hypnorum")
n_bees <- data.frame(n_ind = unlist(n_bees_list), species = names(n_bees_list))
n_labels <- merge(spp_sum, n_bees)

# version for grant BBSRC grant app

m_ratio_spp$species2 <- factor(m_ratio_spp$species, levels = c("B. hypnorum", "B. hortorum", "B. lapidarius", "B. pascuorum", "B. ruderatus", "B. terrestris"))

mrat_plot_app <- ggplot(data = m_ratio_spp, aes(x = species2, y = m))+
                    geom_boxplot()  +
                    geom_boxplot(data = m_ratio_spp[m_ratio_spp$species2 == "B. hypnorum",], aes(x = species2, y = m), fill = "grey")  +
                    #stat_summary(fun.y = mean, geom = "point", shape = 17, size = 2)+
                    ylab("M ratio") +
                    xlab("Species") +
                    geom_text(data = n_labels, mapping = aes(x = species, y = mean - 0.027), label = paste(n_labels$n, paste("(", n_labels$n_ind, ")", sep = ""), sep = " ")) +
                    theme_bw()+
                    theme(panel.grid.major = element_blank())

mrat_plot <- ggplot(data = m_ratio_spp, aes(x = species, y = m))+
  geom_boxplot()  +
  stat_summary(fun.y = mean, geom = "point", shape = 17, size = 2)+
  ylab("M ratio") +
  xlab("Species") +
  geom_text(data = n_labels, mapping = aes(x = species, y = mean - 0.027), label = paste(n_labels$n, paste("(", n_labels$n_ind, ")", sep = ""), sep = " ")) +
  theme_bw()+
  theme(panel.grid.major = element_blank())
# one way aov
oneway.test(m ~ species, data = m_ratio_spp)
mrat_plot
##prep new object 4 bottleneck
worker_gtypes <- read.csv("./data out/g.worker.mat_all_4colony.csv")
worker_gtypes$sample <- as.character(worker_gtypes$sample)
worker_gtypes <- worker_gtypes[worker_gtypes$sample %in% workers,]
write.csv(worker_gtypes, "./data out/filtered.worker.gtypes.csv", row.names = FALSE)

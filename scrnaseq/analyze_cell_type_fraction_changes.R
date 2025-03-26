#### Cell type fraction analysis in snRNA-seq ####
##Adapted from Lai et al. 2024 Cell Rep and https://github.com/JHarrisonEcoEvo/DMM/blob/master/reanalysisLungsHMC.R

options(stringsAsFactors = FALSE)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bayestestR)
library(reshape2)
library(gtools)
library(shinystan)
library(dplyr)
library(StanHeaders)

# # 0. Install RStan and bayestestR ---------------------------------------------
# install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
# install.packages("bayestestR")

# 1. Run model ------------------------------------------------------------------
df <- read.table("n_cells.neuronal.txt",sep="\t",header=T)

# Compare ALS/FTD with Control by region ####
df$region_clinical_TDP43 <- factor(df$region_clinical_TDP43,
                                   c("BA6_Control_normal","BA6_ALS_normal","BA6_ALS_TDP43",
                                     "PFC_Control_normal","PFC_FTD_normal","PFC_FTD_TDP43"))
df.ba6 <- df[df$brain_region=="BA6",c("orig.ident","celltype_detailed","n","region_clinical_TDP43")]
df.pfc <- df[df$brain_region=="PFC",c("orig.ident","celltype_detailed","n","region_clinical_TDP43")]
df.control <- df[df$clinical=="Control",c("orig.ident","celltype_detailed","n","region_clinical_TDP43")]
df.ba6 <- dcast(df.ba6, orig.ident + region_clinical_TDP43 ~ celltype_detailed, value.var = "n", fill=0)
df.pfc <- dcast(df.pfc, orig.ident + region_clinical_TDP43 ~ celltype_detailed, value.var = "n", fill=0)
df.control <- dcast(df.control, orig.ident + region_clinical_TDP43 ~ celltype_detailed, value.var = "n", fill=0)

#order by region_clinical_TDP43
df.ba6 <- df.ba6[order(df.ba6$region_clinical_TDP43),]
df.pfc <- df.pfc[order(df.pfc$region_clinical_TDP43),]
df.control <- df.control[order(df.control$region_clinical_TDP43),]

#get grouping indices (each start and end pair is a group), N = number of groups
starts.ba6 <- c(1,3,4)
ends.ba6 <- c(2,3,5)
starts.pfc <- c(1,3,5)
ends.pfc <- c(2,4,6)
starts.control <- c(1,3)
ends.control <- c(2,4)

#compile model
DM <- rstan::stan_model("DMM_HMC/DM.stan", 
                        model_name="DM")
#run model
fitstan_HMC.ba6 <-rstan::sampling(DM,
                                  data=list("datamatrix"=df.ba6[,3:ncol(df.ba6)],
                                            "nsamples"=nrow(df.ba6[,3:ncol(df.ba6)]),
                                            "ncelltypes"=ncol(df.ba6[,3:ncol(df.ba6)]),
                                            "N"=3,
                                            "start" = starts.ba6,
                                            "end" = ends.ba6),
                                  chains=16,
                                  control = list(max_treedepth = 15),
                                  warmup=5000,
                                  iter=10000,
                                  thin=2,
                                  algorithm="NUTS",
                                  cores=8,
                                  seed=123,
                                  pars<-c("pi")
)

fitstan_HMC.pfc <-rstan::sampling(DM,
                                  data=list("datamatrix"=df.pfc[,3:ncol(df.pfc)],
                                            "nsamples"=nrow(df.pfc[,3:ncol(df.pfc)]),
                                            "ncelltypes"=ncol(df.pfc[,3:ncol(df.pfc)]),
                                            "N"=3,
                                            "start" = starts.pfc,
                                            "end" = ends.pfc),
                                  chains=16,
                                  control = list(max_treedepth = 15),
                                  warmup=5000,
                                  iter=10000,
                                  thin=2,
                                  algorithm="NUTS",
                                  cores=8,
                                  seed=123,
                                  pars<-c("pi")
)


# 2. Test model parameters (num_chain, num_iter) --------------------------------
HDIofMCMC <- function(sampleVec) {
  # Computes highest density interval from a sample of representative values,
  # estimated as shortest credible interval.
  # Arguments:
  # sampleVec
  # is a vector of representative values from a probability distribution.
  # credMass
  # is a scalar between 0 and 1, indicating the mass within the credible
  # interval that is to be estimated.
  # Value:
  # HDIlim is a vector containing the limits of the HDI
  hdi <- ci(sampleVec, method = "HDI") # default 95% CI
  HDIlim <- c(hdi$CI_low, hdi$CI_high)
  return(HDIlim)
}
calc_certain_ratios <- function(ratios_HMC, dimension) {
  positives <- vector()
  negatives <- vector()
  for(i in 1:dim(ratios_HMC)[dimension]){
    hdi <- ci(ratios_HMC[,i], method = "HDI") # default 95% CI
    HDIlim <- c(hdi$CI_low, hdi$CI_high)
    if(between(0, hdi$CI_low, hdi$CI_high) == 'FALSE'){
      positives <- c(positives, i)
    }else{
      negatives <- c(negatives, i)
    }
  }
  return(list(positives = positives, negatives = negatives))
}
# function for plotting
plotr <- function(x, y, z, title, x_labels){
  par(mar=c(10, 5, 2, 2), cex.axis=1)
  segs <- apply(x, y, HDIofMCMC) #lower-bound = segs[1,], higher-bound = segs[2,]
  mylim <- ceiling(max(abs(segs)))
  plot(1:notus, apply(x, y, mean), cex = 1, ylim = c(-mylim,mylim),
       ylab = "Log2 ratio in relative abundance", xlab = " ", main = title,
       cex.lab = 1, cex.main=1, pch = 16, col = colorPoints, las = 2, xaxt="n") +
    axis(1, at=1:notus, labels=x_labels, las=2)
  #mtext(text="Cell type", side=1, line=13.5, cex=1.5)
  abline(h = 0, col = "black", lty=2)
  colorLines <- rep("black", notus)
  colorLines[z$positives] <- "red"
  segments(1:notus, segs[1,], 1:notus, segs[2,], col = colorLines)
}

# Define cell type order
ordered_celltypes <- c("L2-3 IT","L4 IT","L5 Bentz VEN","L5 corticobulbar tract","L5 IT",
                       "L6 CT","L6 IT","L6 IT Car3","L6b","Ch IN","LAMP5+ IN","PV+ IN",
                       "SST+ IN","SST+ LAMP5+ IN","VIP+ IN","VIP+like IN")

##BA6
disease <- "ALS"
region <- "BA6"
my_order <- match(ordered_celltypes, colnames(df.ba6)[3:ncol(df.ba6)])
est.pi <- extract(fitstan_HMC.ba6,"pi")
ratios_HMC.TDPneg_TDPpos <- log2(est.pi$pi[,3,]/est.pi$pi[,2,])
# reorder cell type
ratios_HMC.TDPneg_TDPpos <- ratios_HMC.TDPneg_TDPpos[,my_order]
notus <- dim(ratios_HMC.TDPneg_TDPpos)[2] # number of celltypes
colorPoints <- rep("black", notus)
plotr(x = ratios_HMC.TDPneg_TDPpos, y = 2, z = calc_certain_ratios(ratios_HMC.TDPneg_TDPpos, 2),
      title = paste0(disease,"_TDP43- vs. ",disease,"_TDP43+"),
      x_labels = ordered_celltypes)

##PFC
disease <- "FTD"
region <- "PFC"
my_order <- match(ordered_celltypes, colnames(df.pfc)[3:ncol(df.pfc)])
est.pi <- extract(fitstan_HMC.pfc,"pi")
ratios_HMC.TDPneg_TDPpos <- log2(est.pi$pi[,3,]/est.pi$pi[,2,])
# reorder cell type
ratios_HMC.TDPneg_TDPpos <- ratios_HMC.TDPneg_TDPpos[,my_order]
notus <- dim(ratios_HMC.TDPneg_TDPpos)[2] # number of celltypes
colorPoints <- rep("black", notus)
plotr(x = ratios_HMC.TDPneg_TDPpos, y = 2, z = calc_certain_ratios(ratios_HMC.TDPneg_TDPpos, 2),
      title = paste0(disease,"_TDP43- vs. ",disease,"_TDP43+"),
      x_labels = ordered_celltypes)

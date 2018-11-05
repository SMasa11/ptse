#Monte Carlo DGP

library(readstata13)
library(ggplot2)
library(Rcpp)
library(splines)

##Initialization process
loadTrueData <- function()
{
  #data
  tbl <- read.dta13('./tests/testthat/CreponEtAl.dta')
  RCpp <- TRUE

  # Set Seed
  set.seed(23762, kind = NULL, normal.kind = NULL)

  # take log
  tbl$logoutput <- log(tbl$Y+((tbl$Y)^2+1)^(1/2))
  tbl$Yb <- tbl$output_total_bl
  tbl$Yb <- log(tbl$Yb+((tbl$Yb)^2+1)^(1/2))
  tbl$Yb[tbl$Yb <= 0] <- NA

  tbl$self_empl_bl[tbl$self_empl_bl == 0] <- NA
  tbl$samplemodel[tbl$samplemodel == 0] <- NA
  tbl <- na.omit(tbl)

  #tie-breaking
  tiebreak <- runif(length(tbl$logoutput[tbl$logoutput > 0]),-0.001,0.001)
  tbl$logoutput[tbl$logoutput > 0] <- tbl$logoutput[tbl$logoutput > 0] + tiebreak
  tiebreakalt <- runif(length(tbl$Yb[tbl$Yb > 0]),-0.001,0.001)
  tbl$Yb[tbl$Yb > 0] <- tbl$Yb[tbl$Yb > 0] + tiebreakalt

  #D variable
  tbl$client <- tbl$client + 1
  tbl$treatment <- tbl$T

  #cluster
  ## transform into integers
  tbl$demi_paire <- as.numeric(tbl$demi_paire)*10
  tbl$demi_paire <- factor(tbl$demi_paire)
  #tbl$cluster[tbl$cluster == 310] <- NA
  #tbl$cluster[tbl$cluster == 410] <- NA
  tbl <- na.omit(tbl)

  main <- data.frame(logoutput=tbl$logoutput,treatment=tbl$treatment,client=tbl$client,demi_paire=tbl$demi_paire,Yb=tbl$Yb,head_age_bl=tbl$head_age_bl,members_resid_bl=tbl$members_resid_bl,act_livestock_bl = tbl$act_livestock_bl,act_business_bl=tbl$act_business_bl,borrowed_total_bl=tbl$borrowed_total_bl,ccm_resp_activ=tbl$ccm_resp_activ,other_resp_activ=tbl$other_resp_activ,ccm_resp_activ_d=tbl$ccm_resp_activ_d,nadults_resid_bl=tbl$nadults_resid_bl)
  main <- na.omit(main)
  return(main)
}

runEstimateReporting <- function()
{
  options(warn=0)
  #fml <- "logoutput ~ head_age_bl + members_resid_bl | act_livestock_bl + act_business_bl + borrowed_total_bl + ccm_resp_activ + other_resp_activ + ccm_resp_activ_d + nadults_resid_bl"
  fml <- "logoutput ~ head_age_bl + members_resid_bl | act_livestock_bl + act_business_bl + borrowed_total_bl + nadults_resid_bl"
  numKnots <- 1
  numOrder <- 1
  nBoot <- 3

  Ytype <- "Yb"
  main <- loadTrueData()
  reporting <- meanEstimate(plotBeg = 0.25,plotEnd = 0.95,plotBy=0.05,df=main,dDim=2,treatment="treatment",posttreat="client",numKnots=numKnots,numOrder=numOrder,formula=fml,Ytype=Ytype)

  return(reporting)
}

runEstimateCI <- function(rpt)
{
  options(warn=0)
  fml <- "logoutput ~ head_age_bl + members_resid_bl | act_livestock_bl + act_business_bl + borrowed_total_bl + nadults_resid_bl"
  #fml <- "logoutput ~ head_age_bl + members_resid_bl | act_livestock_bl + act_business_bl + borrowed_total_bl + ccm_resp_activ + other_resp_activ + ccm_resp_activ_d + nadults_resid_bl"
  additiveT <- ""
  additiveSplineT <- FALSE
  numKnots <- 1
  numOrder <- 1
  nBoot <- 301

  link <- "logit"
  Ytype <- "Yalt"
  typeSieve <- "spline"
  discY <- FALSE
  main <- loadTrueData()
  CIs <- bootstrapProc(plotBeg = 0.25,plotEnd = 0.95,plotBy=0.05,df=main,dDim=2,treatment="treatment",posttreat="client",meanReporting=rpt,nBoot=nBoot,clusterInference=TRUE,cluster="demi_paire",numKnots=numKnots,numOrder=numOrder,formula=fml,Ytype="Yb", showDisplayMean = TRUE, saveFile = TRUE)
  return(CIs)
}

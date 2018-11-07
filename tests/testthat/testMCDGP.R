#Monte Carlo DGP

library(readstata13)
library(ggplot2)
library(Rcpp)
library(splines)

##Initialization process
generateMonteCarloDGP <- function(nc)
{
  #-----------------------
  # Mock Data Generation
  #----------------------

  #case 1
  n <- 20

  #case 2
  #nc <- 100

  testData <- function(n)
  {
    nu <- runif(n,0,1)
    nu1 <- runif(n,0,1)
    nub <- runif(n,0,1)
    nu0 <- runif(n,0,1)
    eta <- runif(n,0,1)
    xiW1 <- runif(n,0,1)
    xiW2 <- runif(n,0,1)
    xiW3 <- runif(n,0,1)
    clusterShock <- runif(1,0,1)

    U1 <- 0.5*nu + 0.4*nu1 + 0.1*clusterShock
    U0 <- 0.7*nu + 0.2*nu0 + 0.1*clusterShock
    Ub <- 0.7*nu + 0.2*nub + 0.1*clusterShock

    #X under T == 1
    X1 <- (nu*0.8 + eta*0.2 >= 0.5) + 1
    #X under T == 0
    X0 <- 1 #(nu*0.8 + eta*0.2 >= 0.8) + 1

    WDisc1 <- (xiW1 >= 0.5)
    WCont1 <- xiW2
    WCont2 <- xiW3

    Y1 <- 5.1 + 0.2*WDisc1 + 0.8*WCont1 + 0.3*WDisc1*WCont1 + 0.7*(WCont1*WCont2) + 0.5*WCont2 + qnorm(U1)*1.5
    Y0 <- 5.0 + 0.15*WDisc1 + 0.95*WCont1 + 0.12*WDisc1*WCont1 + 0.6*(WCont1*WCont2) + 0.2*WCont2 + qnorm(U0)*0.8
    Yb <- 4.9 + 0.28*WDisc1 + 0.7*WCont1 - 0.2*WDisc1*WCont1 + 0.1*(WCont1*WCont2) + 0.1*WCont2 + qnorm(Ub)*0.5

    WDisc1 <- factor(WDisc1)

    trueData <- data.frame(Y1=Y1,Y0=Y0,Yb=Yb,X1=X1,X0=X0,WD1=WDisc1,WC1=WCont1,WC2=WCont2)
    return(trueData)
  }

  testDataObs <- function(trueData)
  {
    T <- trueData$T
    Y <- T*trueData$Y1 + (1 - T)*trueData$Y0
    Yb <- trueData$Yb
    X <- T*trueData$X1 + (1 - T)*trueData$X0
    WD1 <- trueData$WD1
    WC1 <- trueData$WC1
    WC2 <- trueData$WC2
    main <- data.frame(Y=Y,T=T,WD1=WD1,WC1=WC1,WC2=WC2,X=X,Yb=Yb,cluster=trueData$cluster)

    # print(mean(main$Y[main$T==1]))
    # print(mean(main$Y[main$T==0]))
    # print(mean(main$Y[main$T==1]) - mean(main$Y[main$T==0]))
    # print(mean(main$X[main$T==1]))
    # print(mean(main$T))

    main <- na.omit(main)
    return(main)
  }

  for (c in 1:nc)
  {
    clusterData <- testData(n = n)
    T <- (c %% 2)
    cluster <- c
    clusterData <- cbind(clusterData,T)
    clusterData <- cbind(clusterData,cluster)
    if (c == 1)
    {
      trueDataCluster <- clusterData
    } else {
      trueDataCluster <- rbind(trueDataCluster,clusterData)
    }
  }
  main <- testDataObs(trueData=trueDataCluster)
  return(main)
}


testDataDGP <- function(n)
{
  nu <- runif(n,0,1)
  nu1 <- runif(n,0,1)
  nub <- runif(n,0,1)
  nu0 <- runif(n,0,1)
  eta <- runif(n,0,1)
  xiW1 <- runif(n,0,1)
  xiW2 <- runif(n,0,1)
  xiW3 <- runif(n,0,1)
  clusterShock <- runif(n,0,1)

  U1 <- 0.5*nu + 0.4*nu1 + 0.1*clusterShock
  U0 <- 0.7*nu + 0.2*nu0 + 0.1*clusterShock
  Ub <- 0.7*nu + 0.2*nub + 0.1*clusterShock

  #X under T == 1
  X1 <- (nu*0.8 + eta*0.2 >= 0.5) + 1
  #X under T == 0
  X0 <- 1 #(nu*0.8 + eta*0.2 >= 0.8) + 1

  WDisc1 <- (xiW1 >= 0.5)
  WCont1 <- xiW2
  WCont2 <- xiW3

  Y1 <- 5.1 + 0.2*WDisc1 + 0.8*WCont1 + 0.3*WDisc1*WCont1 + 0.7*(WCont1*WCont2) + 0.5*WCont2 + qnorm(U1)*1.5
  Y0 <- 5.0 + 0.15*WDisc1 + 0.95*WCont1 + 0.12*WDisc1*WCont1 + 0.6*(WCont1*WCont2) + 0.2*WCont2 + qnorm(U0)*0.8
  Yb <- 4.9 + 0.28*WDisc1 + 0.7*WCont1 - 0.2*WDisc1*WCont1 + 0.1*(WCont1*WCont2) + 0.1*WCont2 + qnorm(Ub)*0.5

  WDisc1 <- factor(WDisc1)

  trueData <- data.frame(Y1=Y1,Y0=Y0,Yb=Yb,X1=X1,X0=X0,WD1=WDisc1,WC1=WCont1,WC2=WCont2)
  return(trueData)
}


runEstimateReporting <- function()
{
  # Set Seed
  set.seed(23762, kind = NULL, normal.kind = NULL)
  main <- generateMonteCarloDGP(nc = 25)

  plotBeg <- 0.2
  plotEnd <- 0.8
  plotBy  <- 0.05
  alpha <- 0.05
  perc <- seq(plotBeg,plotEnd,by=plotBy)

  RCpp <- TRUE
  nGridsFine <- 1000
  nGrids <- 30

  options(warn=0)

  additiveT <- "" #"additive"
  additiveSplineT <- FALSE
  numKnots <- 2
  numOrder <- 1

  Ytype <- "Yb"

  main <- na.omit(main)
  fml <- "Y ~ WC1 + WC2 | WD1 "
  reporting <- meanEstimate(df=main,dDim=2,treatment="T",posttreat="X",numKnots=numKnots,formula=fml,numOrder=numOrder,Ytype = Ytype)
  return(reporting)
}

runEstimateCI <- function(rpt)
{
  # Set Seed
  set.seed(23762, kind = NULL, normal.kind = NULL)

  plotBeg <- 0.2
  plotEnd <- 0.8
  plotBy  <- 0.05
  alpha <- 0.05
  perc <- seq(plotBeg,plotEnd,by=plotBy)

  RCpp <- TRUE
  nGridsFine <- 1000
  nGrids <- 30

  options(warn=0)

  additiveT <- "" #"additive"
  additiveSplineT <- FALSE
  numKnots <- 2
  numOrder <- 1
  nBoot <- 3 #301

  link <- "logit"
  Ytype <- "Yb"
  typeSieve <- "spline"
  discY <- FALSE

  main <- generateMonteCarloDGP(nc = 25)
  main <- na.omit(main)
  fml <- "Y ~ WC1 + WC2 | WD1 "
  #rpt <- meanEstimate(df=main,xType=2,plotBeg = 0.25,plotEnd = 0.75,plotBy=0.05,numDisc=numDisc,numCont=numCont,numQuant=numQuant,addDummy="",dumAllAdditive=additiveT,RCpp=RCpp,formula=fml,typeSieve=typeSieve,Ytype=Ytype,discY=discY,numOrder=numOrder,additiveSpline=additiveSplineT,logToLevel=FALSE,link=link)
  vectorCIs <- bootstrapProc(df=main,dDim=2,treatment="T",posttreat="X",meanReporting=rpt,nBoot=nBoot,Ytype=Ytype,clusterInference=TRUE,cluster="cluster",numKnots = numKnots,formula=fml,numOrder=numOrder,showDisplayMean = FALSE,saveFile = FALSE)
  return(vectorCIs)
}


runCombined <- function()
{
  # Set Seed
  set.seed(23762, kind = NULL, normal.kind = NULL)

  plotBeg <- 0.2
  plotEnd <- 0.8
  plotBy  <- 0.05
  alpha <- 0.05
  perc <- seq(plotBeg,plotEnd,by=plotBy)

  RCpp <- TRUE
  nGridsFine <- 1000
  nGrids <- 30

  options(warn=0)

  additiveT <- "" #"additive"
  additiveSplineT <- FALSE
  numKnots <- 2
  numOrder <- 1
  nBoot <- 3 #301

  link <- "logit"
  Ytype <- "Yb"
  typeSieve <- "spline"
  discY <- FALSE

  main <- generateMonteCarloDGP(nc = 25)
  main <- na.omit(main)
  fml <- "Y ~ WC1 + WC2 | WD1 "
  vectorCIs <- ptsEst(df=main,dDim=2,treatment="T",posttreat="X",nBoot=nBoot,Ytype=Ytype,clusterInference=TRUE,cluster="cluster",numKnots = numKnots,formula=fml,numOrder=numOrder,showDisplayMean = FALSE,saveFile = FALSE)
  return(vectorCIs)
}

runCoverage <- function()
{
  # Set Seed
  set.seed(23762, kind = NULL, normal.kind = NULL)

  plotBeg <- 0.25
  plotEnd <- 0.75
  plotBy  <- 0.05
  alpha <- 0.05
  perc <- seq(plotBeg,plotEnd,by=plotBy)

  trueDataCluster <- testDataDGP(n = 100000)

  #for continuation from termination after 187, for NC = 150
  #set.seed(384722, kind = NULL, normal.kind = NULL)

  #for continuation from termination after 225, for NC = 150
  #set.seed(1846841, kind = NULL, normal.kind = NULL)

  #for continuation from termination after 233, for NC = 150
  #set.seed(48645, kind = NULL, normal.kind = NULL)

  #for continuation from termination after 246, for NC = 150
  #set.seed(3874, kind = NULL, normal.kind = NULL)

  #for continuation from termination after 257, for NC = 150
  #set.seed(415687, kind = NULL, normal.kind = NULL)

  #for continuation from termination after 289, for NC = 150
  #set.seed(54865, kind = NULL, normal.kind = NULL)

  #for continuation from termination after 337, for NC = 150
  set.seed(54654, kind = NULL, normal.kind = NULL)

  cdfY11 <- ecdf(trueDataCluster$Y1[trueDataCluster$X1 == 1])
  qY11 <- quantile(cdfY11,perc)
  cdfY01 <- ecdf(trueDataCluster$Y0[trueDataCluster$X1 == 1])
  qY01 <- quantile(cdfY01,perc)
  cdfY12 <- ecdf(trueDataCluster$Y1[trueDataCluster$X1 == 2])
  qY12 <- quantile(cdfY12,perc)
  cdfY02 <- ecdf(trueDataCluster$Y0[trueDataCluster$X1 == 2])
  qY02 <- quantile(cdfY02,perc)

  qdif01 <- qY11 - qY01
  plotData <- data.frame(x=perc,y=qdif01,case=1)
  qdif02 <- qY12 - qY02
  plotData <- rbind(plotData,data.frame(x=perc,y=qdif02,case=2))

  trueCATE1 <- mean(trueDataCluster$Y1[trueDataCluster$X1 == 2]) - mean(trueDataCluster$Y0[trueDataCluster$X1 == 2])
  trueCATE0 <- mean(trueDataCluster$Y1[trueDataCluster$X1 == 1]) - mean(trueDataCluster$Y0[trueDataCluster$X1 == 1])
  trueQTE1 <- qdif02
  trueQTE0 <- qdif01

  CATE1In <- 324
  CATE2In <- 318
  QTE1In <- 317
  QTE2In <- 319
  options(warn=0)

  numKnots <- 1
  numOrder <- 1
  nBoot <- 101 #301

  naCounter <- 0
  for (i in 338:1000)
  {
    main <- generateMonteCarloDGP(nc = 150)
    main <- na.omit(main)
    fml <- "Y ~ WC1 + WC2 | WD1 "
    rpt <- meanEstimate(df=main,dDim=2,treatment="T",posttreat="X",numKnots=numKnots,formula=fml,numOrder=numOrder)
    vectorCIs <- bootstrapProc(df=main,dDim=2,treatment="T",posttreat="X",meanReporting=rpt,nBoot=nBoot,clusterInference=TRUE,cluster="cluster",numKnots = numKnots,formula=fml,numOrder=numOrder,showDisplayMean = FALSE,saveFile = FALSE)
    if (max(is.na(vectorCIs)))
    {
      naCounter <- naCounter + 1
    } else {
      CATECI <- vectorCIs[2:(2+vectorCIs[1] - 1)]
      indexStart <- vectorCIs[1]+3
      QTECI <- vectorCIs[indexStart:length(vectorCIs)]

      print(i)
      cat("CATE0: true", trueCATE0," CI: ",CATECI[1],CATECI[2], "\n")
      if ((trueCATE0 >= CATECI[1]) & (trueCATE0 <= CATECI[2]))
      {
        CATE1In <- CATE1In + 1
      }

      cat("CATE1: true", trueCATE1," CI: ",CATECI[3],CATECI[4], "\n")
      if ((trueCATE1 >= CATECI[3]) & (trueCATE1 <= CATECI[4]))
      {
        CATE2In <- CATE2In + 1
      }

      flag <- 1
      QTE1LB <- QTECI[1:(length(QTECI)/4)]
      QTE1UB <- QTECI[(length(QTECI)/4 + 1):(2*length(QTECI)/4)]
      cat("QTE0: true", trueQTE0, "\n")
      cat("QTE0LB:", QTE1LB, "\n")
      cat("QTE0UB:", QTE1UB, "\n")
      cat("QTE0LB Check:", trueQTE0 - QTE1LB, "\n")
      cat("QTE0UB Check:", QTE1UB - trueQTE0, "\n")
      for (j in seq(1,length(QTECI)/4,by=1))
      {
        if ((trueQTE0[j] < QTE1LB[j]) | (trueQTE0[j] > QTE1UB[j]))
        {
          flag <- 0
        }
      }
      QTE1In <- QTE1In + flag

      QTE2LB <- QTECI[(2*length(QTECI)/4 + 1):(3*length(QTECI)/4)]
      QTE2UB <- QTECI[(3*length(QTECI)/4 + 1):(4*length(QTECI)/4)]
      cat("QTE1: true", trueQTE1, "\n")
      cat("QTE1LB:", QTE2LB, "\n")
      cat("QTE1UB:", QTE2UB, "\n")
      cat("QTE1LB Check:", trueQTE1 - QTE2LB, "\n")
      cat("QTE1UB Check:", QTE2UB - trueQTE1, "\n")
      flag <- 1
      for (j in seq(1,length(QTECI)/4,by=1))
      {
        if ((trueQTE1[j] < QTE2LB[j]) | (trueQTE1[j] > QTE2UB[j]))
        {
          flag <- 0
        }
      }
      QTE2In <- QTE2In + flag
    }
    cat("i,  CATE1,  CATE2,  QTE1,  QTE2 \n")
    cat(i-naCounter,CATE1In," ",CATE2In," ", QTE1In," ",QTE2In,"\n")
    cat(i-naCounter,CATE1In/(i-naCounter)," ",CATE2In/(i-naCounter)," ", QTE1In/(i-naCounter)," ",QTE2In/(i-naCounter),"\n")
    # cat(i," ",CATE1In," ",CATE2In," ",QTE1In," ",QTE2In,"\n")
    # cat(i," ",CATE1In/i," ",CATE2In/i," ",QTE1In/i," ",QTE2In/i,"\n")
  }
  options(warn=1)
}


generateMonteCarloDGPNoDisc <- function(nc)
{
  n <- 20
  testDataNoDisc <- function(n)
  {
  nu <- runif(n,0,1)
  nu1 <- runif(n,0,1)
  nub <- runif(n,0,1)
  nu0 <- runif(n,0,1)
  eta <- runif(n,0,1)
  xiW1 <- runif(n,0,1)
  xiW2 <- runif(n,0,1)
  xiW3 <- runif(n,0,1)
  clusterShock <- runif(1,0,1)

  U1 <- 0.5*nu + 0.4*nu1 + 0.1*clusterShock
  U0 <- 0.7*nu + 0.2*nu0 + 0.1*clusterShock
  Ub <- 0.7*nu + 0.2*nub + 0.1*clusterShock

  #X under T == 1
  X1 <- (nu*0.8 + eta*0.2 >= 0.5) + 1
  #X under T == 0
  X0 <- 1 #(nu*0.8 + eta*0.2 >= 0.8) + 1

  WCont3 <- xiW1
  WCont1 <- xiW2
  WCont2 <- xiW3

  Y1 <- 5.1 + 0.2*WCont3 + 0.8*WCont1 + 0.3*WCont3*WCont1 + 0.7*(WCont1*WCont2) + 0.5*WCont2 + qnorm(U1)*1.5
  Y0 <- 5.0 + 0.15*WCont3 + 0.95*WCont1 + 0.12*WCont3*WCont1 + 0.6*(WCont1*WCont2) + 0.2*WCont2 + qnorm(U0)*0.8
  Yb <- 4.9 + 0.28*WCont3 + 0.7*WCont1 - 0.2*WCont3*WCont1 + 0.1*(WCont1*WCont2) + 0.1*WCont2 + qnorm(Ub)*0.5

  trueData <- data.frame(Y1=Y1,Y0=Y0,Yb=Yb,X1=X1,X0=X0,WC3=WCont3,WC1=WCont1,WC2=WCont2)
  return(trueData)
}

testDataObsNoDisc <- function(trueData)
{
  T <- trueData$T
  Y <- T*trueData$Y1 + (1 - T)*trueData$Y0
  Yb <- trueData$Yb
  X <- T*trueData$X1 + (1 - T)*trueData$X0
  WC3 <- trueData$WC3
  WC1 <- trueData$WC1
  WC2 <- trueData$WC2
  main <- data.frame(Y=Y,T=T,WC3=WC3,WC1=WC1,WC2=WC2,X=X,Yb=Yb,cluster=trueData$cluster)

  # print(mean(main$Y[main$T==1]))
  # print(mean(main$Y[main$T==0]))
  # print(mean(main$Y[main$T==1]) - mean(main$Y[main$T==0]))
  # print(mean(main$X[main$T==1]))
  # print(mean(main$T))

  main <- na.omit(main)
  return(main)
}

for (c in 1:nc)
{
  clusterData <- testDataNoDisc(n = n)
  T <- (c %% 2)
  cluster <- c
  clusterData <- cbind(clusterData,T)
  clusterData <- cbind(clusterData,cluster)
  if (c == 1)
  {
    trueDataCluster <- clusterData
  } else {
    trueDataCluster <- rbind(trueDataCluster,clusterData)
  }
}
main <- testDataObsNoDisc(trueData=trueDataCluster)
return(main)
}

runEstimateReportingNoDisc <- function()
{
  # Set Seed
  set.seed(23762, kind = NULL, normal.kind = NULL)
  main <- generateMonteCarloDGPNoDisc(nc = 25)

  plotBeg <- 0.2
  plotEnd <- 0.8
  plotBy  <- 0.05
  alpha <- 0.05
  perc <- seq(plotBeg,plotEnd,by=plotBy)

  RCpp <- TRUE
  nGridsFine <- 1000
  nGrids <- 30

  options(warn=0)

  additiveT <- "" #"additive"
  additiveSplineT <- FALSE
  numKnots <- 2
  numOrder <- 1

  Ytype <- "Yb"

  main <- na.omit(main)
  fml <- "Y ~ WC1 + WC2 + WC3 | 0"
  reporting <- meanEstimate(df=main,dDim=2,treatment="T",posttreat="X",numKnots=numKnots,formula=fml,numOrder=numOrder,Ytype = Ytype)
  return(reporting)
}

generateMonteCarloDGPNoCont <- function(nc)
{
  n <- 20
  testDataNoCont <- function(n)
  {
    nu <- runif(n,0,1)
    nu1 <- runif(n,0,1)
    nub <- runif(n,0,1)
    nu0 <- runif(n,0,1)
    eta <- runif(n,0,1)
    xiW1 <- runif(n,0,1)
    xiW2 <- runif(n,0,1)
    xiW3 <- runif(n,0,1)
    clusterShock <- runif(1,0,1)

    U1 <- 0.5*nu + 0.4*nu1 + 0.1*clusterShock
    U0 <- 0.7*nu + 0.2*nu0 + 0.1*clusterShock
    Ub <- 0.7*nu + 0.2*nub + 0.1*clusterShock

    #X under T == 1
    X1 <- (nu*0.8 + eta*0.2 >= 0.5) + 1
    #X under T == 0
    X0 <- 1 #(nu*0.8 + eta*0.2 >= 0.8) + 1

    WCont3 <- (xiW1 > 0.5)
    WCont1 <- (xiW2 > 0.5)
    WCont2 <- (xiW3 > 0.5)

    Y1 <- 5.1 + 0.2*WCont3 + 0.8*WCont1 + 0.3*WCont3*WCont1 + 0.7*(WCont1*WCont2) + 0.5*WCont2 + qnorm(U1)*1.5
    Y0 <- 5.0 + 0.15*WCont3 + 0.95*WCont1 + 0.12*WCont3*WCont1 + 0.6*(WCont1*WCont2) + 0.2*WCont2 + qnorm(U0)*0.8
    Yb <- 4.9 + 0.28*WCont3 + 0.7*WCont1 - 0.2*WCont3*WCont1 + 0.1*(WCont1*WCont2) + 0.1*WCont2 + qnorm(Ub)*0.5

    trueData <- data.frame(Y1=Y1,Y0=Y0,Yb=Yb,X1=X1,X0=X0,WC3=WCont3,WC1=WCont1,WC2=WCont2)
    return(trueData)
  }

  testDataObsNoCont <- function(trueData)
  {
    T <- trueData$T
    Y <- T*trueData$Y1 + (1 - T)*trueData$Y0
    Yb <- trueData$Yb
    X <- T*trueData$X1 + (1 - T)*trueData$X0
    WC3 <- trueData$WC3
    WC1 <- trueData$WC1
    WC2 <- trueData$WC2
    main <- data.frame(Y=Y,T=T,WC3=WC3,WC1=WC1,WC2=WC2,X=X,Yb=Yb,cluster=trueData$cluster)

    # print(mean(main$Y[main$T==1]))
    # print(mean(main$Y[main$T==0]))
    # print(mean(main$Y[main$T==1]) - mean(main$Y[main$T==0]))
    # print(mean(main$X[main$T==1]))
    # print(mean(main$T))

    main <- na.omit(main)
    return(main)
  }

  for (c in 1:nc)
  {
    clusterData <- testDataNoCont(n = n)
    T <- (c %% 2)
    cluster <- c
    clusterData <- cbind(clusterData,T)
    clusterData <- cbind(clusterData,cluster)
    if (c == 1)
    {
      trueDataCluster <- clusterData
    } else {
      trueDataCluster <- rbind(trueDataCluster,clusterData)
    }
  }
  main <- testDataObsNoCont(trueData=trueDataCluster)
  return(main)
}


runEstimateReportingNoCont <- function()
{
  # Set Seed
  set.seed(23762, kind = NULL, normal.kind = NULL)
  main <- generateMonteCarloDGPNoCont(nc = 25)

  plotBeg <- 0.2
  plotEnd <- 0.8
  plotBy  <- 0.05
  alpha <- 0.05
  perc <- seq(plotBeg,plotEnd,by=plotBy)

  RCpp <- TRUE
  nGridsFine <- 1000
  nGrids <- 30

  options(warn=0)

  additiveT <- "" #"additive"
  additiveSplineT <- FALSE
  numKnots <- 2
  numOrder <- 1

  Ytype <- "Yb"

  main <- na.omit(main)
  fml <- "Y ~ 0 | WC1 + WC2 + WC3 "
  reporting <- meanEstimate(df=main,dDim=2,treatment="T",posttreat="X",numKnots=numKnots,formula=fml,numOrder=numOrder,Ytype = Ytype)
  return(reporting)
}


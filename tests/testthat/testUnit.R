library(ptse)
library(testthat)
context("unit tests")
test_that("Argument Type Checks ", {
  dataMat <- cbind(runif(100,0,1),runif(100,0,1))
  formula <- "Y ~ W1 + W2 | W3"
  treatment <- "treat"
  posttreat <- "pt"
  typeSieve <- "spline"
  numOrder <- 1
  numKnots <- 1
  dumAllAdditive <- FALSE
  additiveSpline <- FALSE
  nGrids <- 30
  nGridsFine <- 1000
  plotBeg <- 0.25
  plotEnd <- 0.75
  plotBy <- 0.05
  link <- "logit"
  optMute <- TRUE
  Ytype <- "Yb"
  dDim <- 2
  discY <- FALSE
  clusterInference <- FALSE
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY),
               "Argument df must be a data.frame. \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,clusterInference = clusterInference,discY=discY,alpha=0.05),
               "Argument df must be a data.frame. \n",fixed=TRUE)

  dataMat <- data.frame(Y = runif(100,0,1))
  treatment <- 1
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY),
               "Argument treatment must be a single string \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),
               "Argument treatment must be a single string \n",fixed=TRUE)

  treatment <- "treat"
  dDim <- "dummy"
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY),
               "Argument dDim must be an integer. \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),
               "Argument dDim must be an integer. \n",fixed=TRUE)
  dDim <- 2
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha="boo"),"Argument alpha must be a double \n",fixed=TRUE)

  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting="",treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"Argument meanReporting must be a list, generated from meanEstimates",fixed=TRUE)
})

test_that("Argument Value Checks ", {
  dataMat <- data.frame(Y = runif(100,0,1))
  dDim <- 2
  formula <- "Y ~ W1 + W2 | W3"
  treatment <- "treat"
  posttreat <- "pt"
  typeSieve <- "spline"
  numOrder <- 0
  numKnots <- 1
  dumAllAdditive <- FALSE
  additiveSpline <- FALSE
  nGrids <- 30
  nGridsFine <- 1000
  plotBeg <- 0.25
  plotEnd <- 0.75
  plotBy <- 0.05
  link <- "logit"
  optMute <- TRUE
  Ytype <- "Yb"
  discY <- FALSE
  clusterInference <- FALSE
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY),"numOrder must be positive",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"numOrder must be positive",fixed=TRUE)
  numOrder <- 1
  nGridsFine <- -1
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY),"nGridsFine must be positive",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"nGridsFine must be positive",fixed=TRUE)

  nGridsFine <- 1000
  plotBeg <- -1
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY),"plotBeg must be in (0,1) \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"plotBeg must be in (0,1) \n",fixed=TRUE)
  plotEnd <- 1
  plotBeg <- 0.25
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY),"plotEnd must be in (0,1) \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"plotEnd must be in (0,1) \n",fixed=TRUE)

  plotEnd <- 0.25
  plotBeg <- 0.75
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY),"plotBeg must be strictly less than plotEnd \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"plotBeg must be strictly less than plotEnd \n",fixed=TRUE)

  plotEnd <- 0.75
  plotBeg <- 0.25
  plotBy <- 0.9
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY),"plotBy must be positive and less than the range of plots \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"plotBy must be positive and less than the range of plots \n",fixed=TRUE)
  plotBy <- 0.05
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=1.2),"alpha must be in (0,1) \n",fixed=TRUE)

  posttreat <- "DVar"
  treatment <- "T"

  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,disc=discY),"variable specified for posttreat does not exist \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"variable specified for posttreat does not exist \n",fixed=TRUE)
  dataMat <- data.frame(Y = runif(100,0,1), DVar = (runif(100,0,1) > 0.5) + (runif(100,0,1) > 0.5))
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,disc=discY),"variable specified for treatment does not exist \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"variable specified for treatment does not exist \n",fixed=TRUE)
  dataMat <- data.frame(boo = runif(100,0,1), DVar = (runif(100,0,1) > 0.5) + (runif(100,0,1) > 0.5), T = (runif(100,0,1) > 0.5) )
  expect_error(meanEstimate(df=dataMat,dDim=3,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,disc=discY),"variable specified for outcome measure does not exist",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=3,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"variable specified for outcome measure does not exist",fixed=TRUE)

  dataMat <- data.frame(Y = runif(100,0,1), DVar = (runif(100,0,1) > 0.5) + (runif(100,0,1) > 0.5), T = (runif(100,0,1) > 0.5) )
  expect_error(meanEstimate(df=dataMat,dDim=3,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,disc=discY),"variable Yb does not exist in the data.frame",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=3,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"variable Yb does not exist in the data.frame",fixed=TRUE)


  dataMat <- data.frame(Y = runif(100,0,1), DVar = (runif(100,0,1) > 0.5) + (runif(100,0,1) > 0.5), T = (runif(100,0,1) > 0.5), Yb = runif(100,0,1) )
  expect_error(meanEstimate(df=dataMat,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                            numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                            additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,disc=discY),"Number of levels of posttreat does not match with dDim \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=dDim,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"Number of levels of posttreat does not match with dDim \n",fixed=TRUE)
  dataMat <- data.frame(Y = runif(100,0,1), DVar = (runif(100,0,1) > 0.5) + 3*(runif(100,0,1) > 0.5), T = (runif(100,0,1) > 0.5), Yb = runif(100,0,1) )
  expect_error(meanEstimate(df=dataMat,dDim=4,treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                              numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                              additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                              plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,disc=discY),"posttreat: Failure in conversion of factor levels, specify postreat as factor of levels 1,2,...,dDim \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=4,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"posttreat: Failure in conversion of factor levels, specify postreat as factor of levels 1,2,...,dDim \n",fixed=TRUE)

  clusterInference <- TRUE
  expect_error(bootstrapProc(df=dataMat,dDim=4,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"posttreat: Failure in conversion of factor levels, specify postreat as factor of levels 1,2,...,dDim \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=4,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,cluster = "cluster",alpha=0.05),"posttreat: Failure in conversion of factor levels, specify postreat as factor of levels 1,2,...,dDim \n",fixed=TRUE)

  dataMat <- data.frame(Y = runif(100,0,1), DVar = (runif(100,0,1) > 0.5) + 3*(runif(100,0,1) > 0.5), T = (runif(100,0,1) > 0.5), Yb = runif(100,0,1), cluster = factor(round(runif(100,0,1)*30)) )
  clusterInference <- TRUE
  typeSieve <- "Wavelets"
  expect_error(bootstrapProc(df=dataMat,dDim=4,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,alpha=0.05),"posttreat: Failure in conversion of factor levels, specify postreat as factor of levels 1,2,...,dDim \n",fixed=TRUE)
  expect_error(bootstrapProc(df=dataMat,dDim=4,meanReporting=list(1,2),treatment=treatment,posttreat=posttreat,typeSieve =typeSieve,formula=formula,
                             numOrder=numOrder,numKnots=numKnots,dumAllAdditive = dumAllAdditive,
                             additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                             plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,Ytype=Ytype,discY=discY,clusterInference = clusterInference,cluster = "cluster",alpha=0.05),"posttreat: Failure in conversion of factor levels, specify postreat as factor of levels 1,2,...,dDim \n",fixed=TRUE)
}
)

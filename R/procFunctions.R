#' Main Caller for estimating means and bootstrap CIs
#'
#' This function calls meanEstimates() and bootstrapProc() inside.
#'
#' @param df a data.frame, the exact dataset used for the study.
#' @param dDim an integer, positive, the number of values D would take
#' @param treatment a string, variable name for the treatment variable.
#' @param posttreat a string, variable name for the post-treatment variable.
#' @param typeSieve a string, either "formula" or "spline". Default is "spline" taking B-spline for the continuous variables specified in formula
#' @param formula a formula object, for "spline", specify as "outcome ~ continuous | discrete". For "polynomial", use as in other packages like "lm". See the examples for detailed uses.
#' @param numOrder (for typeSieve = "spline") an integer, positive, the maximum order of the spline polynomial
#' @param numKnots (for typeSieve = "spline") an integer, positive, the number of knots for the spline
#' @param dumAllAdditive (for typeSieve = "spline") a boolean, if TRUE, then the dummies are additively included. Default is FALSE and splines are allowed different across discrete variables term.
#' @param additiveSpline (for typeSieve = "spline") a boolean, if TRUE, then the splines are additively combined. Default is FALSE and uses the kronecker product of uni-dimensional splines.
#' @param nGrids an integer, positive, a number of grids for the outcome in the estimation for the conditional cdfs. Default is set to 30. This number crucially affects the computation time as it requires this many of glm estimations separately.
#' @param nGridsFine an integer, positive, a number of grids used for the numerical integration. Default is set to 1000.
#' @param plotBeg a float, a value in (0,1), the lower range of the quantile difference estimate.
#' @param plotEnd a float, a value in (0,1), the upper range of the quantile difference estimate.
#' @param plotBy a float, a positive value, the increment for the quantile difference esimate.
#' @param link a string, a link function used as a binomial family link for glm. Refer to the specification of the link options for the R package glm. Defalt is "logit"
#' @param optMute a boolean, option to toggle off the outputs used during the development. Default is TRUE.
#' @param Ytype a string, either "Yb" or "Y1". Default is "Yb" which uses a proxy variable for Y0. Specification with "Y1" uses Y1 as a proxy for Y0.
#' @param discY a boolean, option to toggle on an experimental feature for the partial identification with discrete proxy variable. Default is FALSE, which assumes absolutely continuous Yb.
#' @param showDisplayMean, a boolean, if TRUE, then it shows a table of CIs intermediary before the whole bootstrap ends. Default is TRUE.
#' @param saveFile, a boolean, if TRUE, it returns the eps files of plots. Default is TRUE.
#' @param clusterInference, a boolean, if TRUE then cluster resampling is used instead of individual resampling
#' @param cluster a string, variable name representing the cluster. Default is "".
#' @param nBoot an integer, positive, a number of bootstrap iterations for the confidence interval construction, plus the true mean estimate. Default is 301.
#' @param alpha a float, a value in (0,1), the size of the test. The default is 0.05.

ptsEst <- function(df,dDim,treatment,posttreat,typeSieve = "spline",formula,numOrder,numKnots,dumAllAdditive = FALSE,additiveSpline = FALSE,nGrids = 30,nGridsFine = 1000,
                   plotBeg = 0.25,plotEnd = 0.75,plotBy=0.05,link="logit",optMute = TRUE,Ytype="Yb",discY=FALSE,showDisplayMean = TRUE,saveFile = TRUE,clusterInference,cluster = "",nBoot=301,alpha=0.05)
{
  report <- meanEstimate(df=df,dDim=dDim,treatment=treatment,posttreat=posttreat,typeSieve=typeSieve,formula=formula,numOrder=numOrder,numKnots=numKnots,dumAllAdditive=dumAllAdditive,
                         additiveSpline = additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
               plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute = optMute,Ytype=Ytype,discY=discY)
  CIs <- bootstrapProc(df=df,dDim=dDim,meanReporting=report,treatment=treatment,posttreat=posttreat,typeSieve = typeSieve,formula=formula,
                        numOrder=numOrder,numKnots=numKnots,dumAllAdditive=dumAllAdditive,additiveSpline=additiveSpline,nGrids = nGrids,nGridsFine = nGridsFine,
                            plotBeg = plotBeg,plotEnd = plotEnd,plotBy=plotBy,link=link,optMute =optMute,clusterInference=clusterInference,cluster =cluster,
                            nBoot=nBoot,alpha=alpha,Ytype=Ytype,discY=discY,showDisplayMean =showDisplayMean,saveFile = saveFile)
  return(CIs)
}


#' Estimation of mean effects and quantile differences
#'
#' This function returns the mean estimates of the average effects and the quantile effects.
#'
#' @param df a data.frame, the exact dataset used for the study.
#' @param dDim an integer, positive, the number of values D would take
#' @param treatment a string, variable name for the treatment variable.
#' @param posttreat a string, variable name for the post-treatment variable.
#' @param typeSieve a string, either "formula" or "spline". Default is "spline" taking B-spline for the continuous variables specified in formula
#' @param formula a formula object, for "spline", specify as "outcome ~ continuous | discrete". For "polynomial", use as in other packages like "lm". See the examples for detailed uses.
#' @param numOrder (for typeSieve = "spline") an integer, positive, the maximum order of the spline polynomial
#' @param numKnots (for typeSieve = "spline") an integer, positive, the number of knots for the spline
#' @param dumAllAdditive (for typeSieve = "spline") a boolean, if TRUE, then the dummies are additively included. Default is FALSE and splines are allowed different across discrete variables term.
#' @param additiveSpline (for typeSieve = "spline") a boolean, if TRUE, then the splines are additively combined. Default is FALSE and uses the kronecker product of uni-dimensional splines.
#' @param nGrids an integer, positive, a number of grids for the outcome in the estimation for the conditional cdfs. Default is set to 30. This number crucially affects the computation time as it requires this many of glm estimations separately.
#' @param nGridsFine an integer, positive, a number of grids used for the numerical integration. Default is set to 1000.
#' @param plotBeg a float, a value in (0,1), the lower range of the quantile difference estimate.
#' @param plotEnd a float, a value in (0,1), the upper range of the quantile difference estimate.
#' @param plotBy a float, a positive value, the increment for the quantile difference esimate.
#' @param link a string, a link function used as a binomial family link for glm. Refer to the specification of the link options for the R package glm. Defalt is "logit"
#' @param optMute a boolean, option to toggle off the outputs used during the development. Default is TRUE.
#' @param Ytype a string, either "Yb" or "Y1". Default is "Yb" which uses a proxy variable for Y0. Specification with "Y1" uses Y1 as a proxy for Y0.
#' @param discY a boolean, option to toggle on an experimental feature for the partial identification with discrete proxy variable. Default is FALSE, which assumes absolutely continuous Yb.
#'
#' @return returns a list, reporting, containing a vector of average effects $reportingV, and a vector of quantile differences $reportingQ.
#'
#' @section Note:
#' typeSieve = "formula" uses the argument formula as a specification of the linear sieve.
#'
#' @usage
#' fml <- "logoutput ~ head_age_bl + members_resid_bl | act_livestock_bl + act_business_bl + borrowed_total_bl + ccm_resp_activ + other_resp_activ + ccm_resp_activ_d + nadults_resid_bl"
#' meanEstimate(df=data,dDim=2,treatment="treatment",posttreat="client",numKnots=1,numOrder=1,formula=fml,Ytype="Yb")
#' # head_age_bl and members_resid_bl are converted to 2-dim B-spline, and the rest will be factorized.
#'
#' fml <- "logoutput ~ head_age_bl + members_resid_bl + head_age_bl*members_resid_bl + act_livestock_bl + act_business_bl + borrowed_total_bl + ccm_resp_activ + other_resp_activ + ccm_resp_activ_d + nadults_resid_bl"
#' meanEstimate(df=data,dDim=2,treatment="treatment",posttreat="client",formula=fml,Ytype="Yb",typeSieve="polynomial")
#' # regard the formula give as the exact polynomial expression to be used for the estimations.
#'
meanEstimate <- function(df,dDim,treatment,posttreat,typeSieve = "spline",formula,numOrder,numKnots,dumAllAdditive = FALSE,additiveSpline = FALSE,nGrids = 30,nGridsFine = 1000,
                         plotBeg = 0.25,plotEnd = 0.75,plotBy=0.05,link="logit",optMute = TRUE,Ytype="Yb",discY=FALSE)
{
  isInteger <- function(int) if (is.numeric(int)) !(int - round(int)) else FALSE

  # argument type checks
  if (!is.data.frame(df)) stop("Argument df must be a data.frame. \n")
  if (!isInteger(dDim)) stop("Argument dDim must be an integer. \n")
  if ((!is.character(treatment)) | (length(treatment) > 1)) stop("Argument treatment must be a single string \n")
  if ((!is.character(posttreat)) | (length(posttreat) > 1)) stop("Argument treatment must be a single string \n")
  if ((!is.character(link)) | (length(link) > 1)) stop("Argument link must be a single string \n")
  if (!isInteger(numOrder)) stop("Argument numOrder must be an integer. \n")
  if (!isInteger(numKnots)) stop("Argument numKnots must be an integer. \n")
  if (!isInteger(nGrids)) stop("Argument nGrids must be an integer. \n")
  if (!isInteger(nGridsFine)) stop("Argument nGridsFine must be an integer. \n")
  if (!((dumAllAdditive == TRUE) | (dumAllAdditive == FALSE))) stop("Argument dumAllAdditive must be a Boolean \n")
  if (!((additiveSpline == TRUE) | (additiveSpline == FALSE))) stop("Argument additiveSpline must be a Boolean \n")
  if (!((optMute == TRUE) | (optMute == FALSE))) stop("Argument optMute must be a Boolean \n")
  if (!((discY == TRUE) | (discY == FALSE))) stop("Argument discY must be a Boolean \n")
  if (!is.double(plotBeg)) stop("Argument plotBeg must be a double \n")
  if (!is.double(plotEnd)) stop("Argument plotEnd must be a double \n")
  if (!is.double(plotBy)) stop("Argument plotBy must be a double \n")
  if (!((Ytype == "Yb") | (Ytype == "Y1"))) stop("Argument Ytype must be either Yb or Y1 \n")

  if (numOrder <= 0) stop("numOrder must be positive \n")
  if (numKnots <= 0) stop("numKnots must be positive \n")
  if (nGrids <= 0) stop("nGrids must be positive \n")
  if (nGridsFine <= 0) stop("nGridsFine must be positive \n")

  if ((plotEnd <= 0) | (plotEnd >= 1)) stop("plotEnd must be in (0,1) \n")
  if ((plotBeg <= 0) | (plotBeg >= 1)) stop("plotBeg must be in (0,1) \n")
  if (plotEnd <= plotBeg) stop("plotBeg must be strictly less than plotEnd \n")
  if ((plotBy <= 0) | plotBy >= (plotEnd - plotBeg)) stop("plotBy must be positive and less than the range of plots \n")

  if (numKnots > 4) stop("optKnots > 4 is not implemented \n")

  if (discY == TRUE) message("Warning: discrete Y is only experimental/not-completed. \n")

  err <- tryCatch(df[,posttreat],error = function(e) e)
  if (is(err,"error")) stop("variable specified for posttreat does not exist \n")

  err <- tryCatch(df[,treatment],error = function(e) e)
  if (is(err,"error")) stop("variable specified for treatment does not exist \n")

  fml <- Formula::as.Formula(formula)
  err <- tryCatch(df[,all.vars(fml)[1]],error = function(e) e)
  if (is(err, "error")) stop("variable specified for outcome measure does not exist")

  if (Ytype == "Yb")
  {
    err <- tryCatch(df[,"Yb"],error = function(e) e)
    if (is(err, "error")) stop("variable Yb does not exist in the data.frame")
  }

  if (!is.factor(df[,posttreat])) df[,posttreat] <- factor(df[,posttreat])

  if (!(length(levels(df[,posttreat])) == dDim))
  {
    stop("Number of levels of posttreat does not match with dDim \n")
  } else {
    if (!(min(levels(df[,posttreat]) == c(1:dDim))))
    {
      if ((min(levels(df[,posttreat]) == c(0:(dDim-1)))) | (min(levels(df[,posttreat]) == c(FALSE,TRUE))))
      {
        levels(df[,posttreat]) <- c(1:dDim)
      } else {
        stop("posttreat: Failure in conversion of factor levels, specify postreat as factor of levels 1,2,...,dDim \n")
      }
    }
  }


  if (!is.factor(df[,treatment])) df[,treatment] <- factor(df[,treatment])

  if (!(length(levels(df[,treatment])) == 2))
  {
    stop("Number of levels of treatment must be 2 \n")
  } else {
    if (!(min(levels(df[,treatment]) == c(0:1))))
    {
      if ((min(levels(df[,treatment]) == c(1:2))) | (min(levels(df[,treatment]) == c(FALSE,TRUE))))
      {
        levels(df[,treatment]) <- c(0:1)
      } else {
        stop("treatment: Failure in conversion of factor levels, specify postreat as factor of levels 1,2,...,dDim \n")
      }
    }
  }

  if (!((typeSieve == "spline") | (typeSieve == "polynomial"))) stop("Argument typeSieve must be either spline or polynomial \n")


  if (dumAllAdditive)
  {
    dumAllAdditive <- "additive"
  } else {
    dumAllAdditive <- ""
  }
  optset <- list("typeSieve" = typeSieve,"optMute" = optMute,"nGrids" = nGrids,"nGridsFine" = nGridsFine,
                 "plotBeg" = plotBeg,"plotEnd" = plotEnd,"plotBy"=plotBy,
                 "dDim"=dDim,"numKnots"=numKnots,
                 "dumAllAdditive"=dumAllAdditive,"discY" = discY,"numOrder"=numOrder,"additiveSpline"=additiveSpline)

  df <- reformatFormula(fml = formula,dataOrigin = df,treatment=treatment,posttreat=posttreat)
  optset["numCont"] <- length(attr(terms(fml,rhs=1),"term.labels"))
  optset["numDisc"] <- length(attr(terms(fml,rhs=2),"term.labels"))
  # mean estimates
  reporting <- mainProcedure(main=df,optset=optset,fml=formula,Ytype=Ytype,link=link)
  return(reporting)
}

#' Bootstrap procedure after the mean estimates are given
#'
#' This function runs bootstrap procedure and returns the confidence intervals for the mean estimates and the uniform confidence bands for the quantile differences. It also produces quantile difference plots.
#'
#' @param df a data.frame, the exact dataset used for the study.
#' @param dDim an integer, positive, the number of values D would take.
#' @param treatment a string, variable name for the treatment variable.
#' @param posttreat a string, variable name for the post-treatment variable.
#' @param meanReporting a list, reporting list produced from the meanEstimate function.
#' @param typeSieve a string, either "polynomial" or "spline". Default is "spline" taking B-spline for the continuous variables specified in formula
#' @param formula a formula object, for "spline", specify as "outcome ~ continuous | discrete". For "polynomial", use as in other packages like "lm". See the examples for detailed uses.
#' @param numOrder (for typeSieve = "spline") an integer, positive, the maximum order of the spline polynomial. It must be the same as previously used for meanEstimate.
#' @param numKnots (for typeSieve = "spline") an integer, positive, the number of knots for the spline. It must be the same as previously used for meanEstimate.
#' @param clusterInference, a boolean, if TRUE then cluster resampling is used instead of individual resampling
#' @param cluster a string, variable name representing the cluster. Default is "".
#' @param nBoot an integer, positive, a number of bootstrap iterations for the confidence interval construction, plus the true mean estimate. Default is 301.
#' @param alpha a float, a value in (0,1), the size of the test. The default is 0.05.
#' @param dumAllAdditive (for typeSieve = "spline") a boolean, if TRUE, then the dummies are additively included. Default is FALSE and splines are allowed different across discrete variables term.
#' @param additiveSpline (for typeSieve = "spline") a boolean, if TRUE, then the splines are additively combined. Default is FALSE and uses the kronecker product of uni-dimensional splines.
#' @param nGrids an integer, positive, a number of grids for the outcome in the estimation for the conditional cdfs. Default is set to 30. This number crucially affects the computation time as it requires this many of glm estimations separately.
#' @param nGridsFine an integer, positive, a number of grids used for the numerical integration. Default is set to 1000.
#' @param plotBeg a float, a value in (0,1), the lower range of the quantile difference estimate.
#' @param plotEnd a float, a value in (0,1), the upper range of the quantile difference estimate.
#' @param plotBy a float, a positive value, the increment for the quantile difference esimate.
#' @param link a string, a link function used as a binomial family link for glm. Refer to the specification of the link options for the R package glm. Defalt is "logit"
#' @param optMute a boolean, option to toggle off the outputs used during the development. Default is TRUE.
#' @param Ytype a string, either "Yb" or "Y1". Default is "Yb" which uses a proxy variable for Y0. Specification with "Y1" uses Y1 as a proxy for Y0.
#' @param discY a boolean, option to toggle on an experimental feature for the partial identification with discrete proxy variable. Default is FALSE, which assumes absolutely continuous Yb.
#' @param showDisplayMean, a boolean, if TRUE, then it shows a table of CIs intermediary before the whole bootstrap ends. Default is TRUE.
#' @param saveFile, a boolean, if TRUE, it returns the eps files of plots. Default is TRUE.
#'
#' @return returns the vector of CIs. The first element is the length of mean effect CIs and the mean CIs, then followed by the length of quantile CIs and the quantile CIs.
#'

bootstrapProc <- function(df,dDim,meanReporting,treatment,posttreat,typeSieve = "spline",formula,numOrder,numKnots,dumAllAdditive=FALSE,additiveSpline=FALSE,nGrids = 30,nGridsFine = 1000,
                          plotBeg = 0.25,plotEnd = 0.75,plotBy=0.05,link="logit",optMute = TRUE,clusterInference,cluster = "",
                          nBoot=301,alpha=0.05,Ytype="Yb",discY=FALSE,showDisplayMean = TRUE,saveFile = TRUE)
{
  isInteger <- function(int) if (is.numeric(int)) !(int - round(int)) else FALSE

  # argument type checks
  if (!is.data.frame(df)) stop("Argument df must be a data.frame. \n")
  if (!is.list(meanReporting)) stop("Argument meanReporting must be a list, generated from meanEstimates")
  if (!isInteger(dDim)) stop("Argument dDim must be an integer. \n")
  if ((!is.character(treatment)) | (length(treatment) > 1)) stop("Argument treatment must be a single string \n")
  if ((!is.character(posttreat)) | (length(posttreat) > 1)) stop("Argument treatment must be a single string \n")
  if ((!is.character(link)) | (length(link) > 1)) stop("Argument link must be a single string \n")
  if (!isInteger(numOrder)) stop("Argument numOrder must be an integer. \n")
  if (!isInteger(numKnots)) stop("Argument numKnots must be an integer. \n")
  if (!isInteger(nGrids)) stop("Argument nGrids must be an integer. \n")
  if (!isInteger(nGridsFine)) stop("Argument nGridsFine must be an integer. \n")
  if (!((dumAllAdditive == TRUE) | (dumAllAdditive == FALSE))) stop("Argument dumAllAdditive must be a Boolean \n")
  if (!((additiveSpline == TRUE) | (additiveSpline == FALSE))) stop("Argument additiveSpline must be a Boolean \n")
  if (!((optMute == TRUE) | (optMute == FALSE))) stop("Argument optMute must be a Boolean \n")
  if (!((discY == TRUE) | (discY == FALSE))) stop("Argument discY must be a Boolean \n")
  if (!((clusterInference == TRUE) | (clusterInference == FALSE))) stop("Argument clusterInference must be a Boolean \n")
  if (!is.double(plotBeg)) stop("Argument plotBeg must be a double \n")
  if (!is.double(plotEnd)) stop("Argument plotEnd must be a double \n")
  if (!is.double(plotBy)) stop("Argument plotBy must be a double \n")
  if (!is.double(alpha)) stop("Argument alpha must be a double \n")
  if (!((Ytype == "Yb") | (Ytype == "Y1"))) stop("Argument Ytype must be either Yb or Y1")

  if (numOrder <= 0) stop("numOrder must be positive")
  if (numKnots <= 0) stop("numKnots must be positive")
  if (nGrids <= 0) stop("nGrids must be positive")
  if (nGridsFine <= 0) stop("nGridsFine must be positive")


  if ((plotEnd <= 0) | (plotEnd >= 1)) stop("plotEnd must be in (0,1) \n")
  if ((plotBeg <= 0) | (plotBeg >= 1)) stop("plotBeg must be in (0,1) \n")
  if (plotEnd <= plotBeg) stop("plotBeg must be strictly less than plotEnd \n")
  if ((plotBy <= 0) | plotBy >= (plotEnd - plotBeg)) stop("plotBy must be positive and less than the range of plots \n")

  if ((alpha <= 0) | (alpha >= 1)) stop("alpha must be in (0,1) \n")

  err <- tryCatch(df[,posttreat],error = function(e) e)
  if (is(err,"error")) stop("variable specified for posttreat does not exist \n")

  err <- tryCatch(df[,treatment],error = function(e) e)
  if (is(err,"error")) stop("variable specified for treatment does not exist \n")


  fml <- Formula::as.Formula(formula)
  err <- tryCatch(df[,all.vars(fml)[1]],error = function(e) e)
  if (is(err, "error")) stop("variable specified for outcome measure does not exist")

  if (Ytype == "Yb")
  {
    err <- tryCatch(df[,"Yb"],error = function(e) e)
    if (is(err, "error")) stop("variable Yb does not exist in the data.frame")
  }
  if (!is.factor(df[,posttreat])) df[,posttreat] <- factor(df[,posttreat])

  if (!(length(levels(df[,posttreat])) == dDim))
  {
    stop("Number of levels of posttreat does not match with dDim \n")
  } else {
    if (!(min(levels(df[,posttreat]) == c(1:dDim))))
    {
      if ((min(levels(df[,posttreat]) == c(0:(dDim-1)))) | (min(levels(df[,posttreat]) == c(FALSE,TRUE))))
      {
        levels(df[,posttreat]) <- c(1:dDim)
      } else {
        stop("posttreat: Failure in conversion of factor levels, specify postreat as factor of levels 1,2,...,dDim \n")
      }
    }
  }


  if (!is.factor(df[,treatment])) df[,treatment] <- factor(df[,treatment])

  if (!(length(levels(df[,treatment])) == 2))
  {
    stop("Number of levels of treatment must be 2 \n")
  } else {
    if (!(min(levels(df[,treatment]) == c(0:1))))
    {
      if ((min(levels(df[,treatment]) == c(1:2))) | (min(levels(df[,treatment]) == c(FALSE,TRUE))))
      {
        levels(df[,treatment]) <- c(0:1)
      } else {
        stop("treatment: Failure in conversion of factor levels, specify postreat as factor of levels 1,2,...,dDim \n")
      }
    }
  }

  if (clusterInference == TRUE)
  {
    err <- tryCatch(df[,cluster],error = function(e) e)
    if (is(err,"error")) stop("variable specified for cluster does not exist \n")
  }

  if (!((typeSieve == "spline") | (typeSieve == "polynomial"))) stop("Argument typeSieve must be either spline or polynomial \n")


  if (dumAllAdditive)
  {
    dumAllAdditive <- "additive"
  } else {
    dumAllAdditive <- ""
  }
  optset <- list("typeSieve" = typeSieve,"optMute" = optMute,"nGrids" = nGrids,"nGridsFine" = nGridsFine,
                 "plotBeg" = plotBeg,"plotEnd" = plotEnd,"plotBy"=plotBy,
                 "dDim"=dDim,"nBoot"=nBoot,"numKnots"=numKnots,
                 "dumAllAdditive" = dumAllAdditive,"discY" = discY,"numOrder"=numOrder,"additiveSpline"=additiveSpline)
  debugFinder <- FALSE
  df <- reformatFormula(fml = formula,dataOrigin = df,treatment=treatment,posttreat=posttreat,cluster=cluster)
  optset["numCont"] <- length(attr(terms(fml,rhs=1),"term.labels"))
  optset["numDisc"] <- length(attr(terms(fml,rhs=2),"term.labels"))

  # extracting means
  if (optset["discY"] == FALSE)
  {
    repV <- cbind(meanReporting["reportingV"][[1]])
  } else {
    repVUB <- cbind(meanReporting["reportingVUB"][[1]])
    repVLB <- cbind(meanReporting["reportingVLB"][[1]])
  }

  # extracting quantiles: concatinate the whole cases for each X value
  if (optset["discY"] == FALSE)
  {
    Q <- meanReporting["reportingQ"][[1]]
    # container for Quantiles
    Qmat <- matrix(unlist(optset["nBoot"])*length(Q),length(Q),unlist(optset["nBoot"]))
    Qmat[,1] <- Q
  } else {
    QUB <- meanReporting["reportingQUB"][[1]]
    QLB <- meanReporting["reportingQLB"][[1]]
    # container for Quantiles
    QUBmat <- matrix(unlist(optset["nBoot"])*length(QUB),length(QUB),unlist(optset["nBoot"]))
    QUBmat[,1] <- QUB
    QLBmat <- matrix(unlist(optset["nBoot"])*length(QLB),length(QLB),unlist(optset["nBoot"]))
    QLBmat[,1] <- QLB
  }
  if (debugFinder)
  {
    print("here")
  }

  # display quantile estimates mean plots
  for (i in 1:unlist(optset["dDim"]))
  {
    if (optset["discY"] == FALSE)
    {
      if (i == 1)
      {
        plotDt <- data.frame(x = seq(unlist(optset["plotBeg"]),unlist(optset["plotEnd"]),by=unlist(optset["plotBy"])),y=Q[1:(length(Q)/unlist(optset["dDim"]))],case=i)
      }
      else
      {
        plotDt <- rbind(plotDt,data.frame(x = seq(unlist(optset["plotBeg"]),unlist(optset["plotEnd"]),by=unlist(optset["plotBy"])),y=Q[(1+length(Q)*(i-1)/unlist(optset["dDim"])):(length(Q)*i/unlist(optset["dDim"]))],case=i))
        print(ggplot(plotDt,aes(x=x,y=y)) + geom_line(aes(color=factor(case))))
      }
    } else {
      if (i == 1)
      {
        plotDt <- data.frame(x = seq(unlist(optset["plotBeg"]),unlist(optset["plotEnd"]),by=unlist(optset["plotBy"])),y=QUB[1:(length(QUB)/unlist(optset["dDim"]))],case=i)
      }
      else
      {
        plotDt <- rbind(plotDt,data.frame(x = seq(unlist(optset["plotBeg"]),unlist(optset["plotEnd"]),by=unlist(optset["plotBy"])),y=QUB[(1+length(QUB)*(i-1)/unlist(optset["dDim"])):(length(QUB)*i/unlist(optset["dDim"]))],case=i))
        print(ggplot(plotDt,aes(x=x,y=y)) + geom_line(aes(color=factor(case))))
      }
    }
  }
  if (debugFinder)
  {
    print("here")
  }

  for (i in 1:unlist(optset["dDim"]))
  {
    if (optset["discY"] == TRUE)
    {
      if (i == 1)
      {
        plotDt <- data.frame(x = seq(unlist(optset["plotBeg"]),unlist(optset["plotEnd"]),by=unlist(optset["plotBy"])),y=QLB[1:(length(QLB)/unlist(optset["dDim"]))],case=i)
      }
      else
      {
        plotDt <- rbind(plotDt,data.frame(x = seq(unlist(optset["plotBeg"]),unlist(optset["plotEnd"]),by=unlist(optset["plotBy"])),y=QLB[(1+length(QLB)*(i-1)/unlist(optset["dDim"])):(length(QLB)*i/unlist(optset["dDim"]))],case=i))
        print(ggplot(plotDt,aes(x=x,y=y)) + geom_line(aes(color=factor(case))))
      }
    }
  }
  if (debugFinder)
  {
    print("here")
  }

  # sample sizes
  n0 <- sum(df$T==0)
  n1 <- sum(df$T==1)
  n <- n0 + n1

  # Cluster bootstrap as iid resampling of clusters with cluster size weight
  ## If it involves cluster, data.frame must contain "cluster" as list of integers
  if (clusterInference == TRUE)
  {

    # get the list of cluster ids
    ## treatment region
    clusterList1 <- unique(sort(df$cluster[df$T==1]))
    ## control region
    clusterList0 <- unique(sort(df$cluster[df$T==0]))

    # fraction for each cluster: container
    #pCluster1 <- c(1:length(clusterList1))*0
    nCluster1 <- c(1:length(clusterList1))*0
    #pCluster0 <- c(1:length(clusterList0))*0
    nCluster0 <- c(1:length(clusterList0))*0
    ## fillin the values
    for (i in 1:length(nCluster1))
    {
      #pCluster1[i] <- sum(df$cluster == clusterList1[i])/sum(df$T == 1)
      nCluster1[i] <- sum(df$cluster == clusterList1[i])
    }
    for (i in 1:length(nCluster0))
    {
      #pCluster0[i] <- sum(df$cluster == clusterList0[i])/sum(df$T == 0)
      nCluster0[i] <- sum(df$cluster == clusterList0[i])
    }
    clList1 <- c(1:length(clusterList1))
    clList0 <- c(1:length(clusterList0))
  }
  # Drawing bootstrapped sample
  ## Draw nBoot
  timeAcc <- 0
  for (i in 1:(unlist(optset["nBoot"])-1))
  {
    startTime <- Sys.time()
    # Resampling procedure
    ## Cluster version
    if (clusterInference == TRUE)
    {
      ### sample from cluster list with
      clDrawn1 <- sample(clList1,size=n1*2,replace = TRUE)
      count <- 1
      sizeBoot1 = nCluster1[clDrawn1[count]]
      ### create new data.frame of the drawn sample
      mainBoot <- df[df$cluster == clusterList1[clDrawn1[count]],]
      #for (looper in 1:(length(clList1)-1))
      while (sizeBoot1 < n1)
      {
        count <- count + 1
        ### add rows of drawn clusters to mainBoot
        mainBoot <- rbind(mainBoot,df[df$cluster == clusterList1[clDrawn1[count]],])
        sizeBoot1 <- sizeBoot1 + nCluster1[clDrawn1[count]]
      }
      #print(sizeBoot1)
      ### continue until drawn sample reaches the original sample size

      ## sampe procedure for T == 0
      clDrawn0 <- sample(clList0,size=n0*2,replace = TRUE)
      count <- 1
      sizeBoot0 = nCluster0[clDrawn0[count]]
      #for (looper in 1:(length(clList0)-1))
      while (sizeBoot0 < n1)
      {
        count <- count + 1
        mainBoot <- rbind(mainBoot,df[df$cluster == clusterList0[clDrawn0[count]],])
        sizeBoot0 <- sizeBoot0 + nCluster0[clDrawn0[count]]
      }
    } else {
      # for non-clustered iid sample
      ## This part needs some modification if original data is generated from stratified sample
      subMain1 <- df[df$T==1,]
      ## sample from T = 1
      mainBoot <- subMain1[sample(nrow(subMain1),size=sum(df$T==1),replace=TRUE),]
      subMain0 <- df[df$T==0,]
      ## sample from T = 0
      mainBoot <- rbind(mainBoot,subMain0[sample(nrow(subMain0),size=sum(df$T==0),replace=TRUE),])
    }
    ## Using the mainBoot data.frame, return the mean estimates
    reporting <- mainProcedure(main=mainBoot,optset=optset,fml=formula,Ytype=Ytype,link=link)

    ## column bind the ith boot for mean effects to repV
    if (optset["discY"] == FALSE)
    {
      repV <- cbind(repV,cbind(reporting["reportingV"][[1]]))
    } else {
      repVUB <- cbind(repVUB,cbind(reporting["reportingVUB"][[1]]))
      repVLB <- cbind(repVLB,cbind(reporting["reportingVLB"][[1]]))
    }
    ## insert quantile estimates to the container
    if (optset["discY"] == FALSE)
    {
      Qmat[,i+1] <- reporting["reportingQ"][[1]]
    } else {
      QUBmat[,i+1] <- reporting["reportingQUB"][[1]]
      QLBmat[,i+1] <- reporting["reportingQLB"][[1]]
    }
    endTime <- Sys.time()
    diffTime <- (endTime - startTime)
    if (attr(diffTime,"units") == "mins")
    {
      timeRound <- as.double(diffTime*60)
    } else {
      timeRound <- as.double(diffTime)
    }
    timeAcc <- timeAcc + timeRound
    if (i == 1)
    {
      message("Expected time ", floor(timeRound*unlist(optset["nBoot"])/60), " mins ", round(timeRound*unlist(optset["nBoot"]) - floor(timeRound*unlist(optset["nBoot"])/60)*60), " secs")
    }
    if ((i %% 50 == 0) | (i == 10))
    {
      message(i," th boot out of ", unlist(optset["nBoot"]))
      if (i == 10)
      {
        timeAvg <- timeAcc/10
      } else {
        timeAvg <- timeAcc/50
      }
      message("Expected remaining time ", floor(timeAvg*(unlist(optset["nBoot"]) - i)/60), " mins ", round(timeAvg*(unlist(optset["nBoot"]) - i) - floor(timeAvg*(unlist(optset["nBoot"]) - i)/60)*60), " secs")
      timeAcc <- 0
      if ((i %% 30 == 0) & (i < unlist(optset["nBoot"])-1))
      {
        if (showDisplayMean)
        {
          if (optset["discY"] == FALSE)
          {
            dump <- displayMeans(repV=repV,optset=optset,n=n,curB=i+1)
          } else {
            dump <- displayMeans(repV=repVUB,optset=optset,n=n,repVLB=repVLB,curB=i+1)
          }
        }
      }
    }
  }
  if (debugFinder)
  {
    print("here")
  }

  # this procedure generate a formatted table for mean, se and CIs for CATE, ATE and ATEPredicted
  if (optset["discY"] == FALSE)
  {
    meanCIs <- displayMeans(repV=repV,optset=optset,n=n,curB=unlist(optset["nBoot"]),showDisplayMean=showDisplayMean)
  } else {
    meanCIs <- displayMeans(repV=repVUB,optset=optset,n=n,repVLB=repVLB,curB=unlist(optset["nBoot"]),showDisplayMean=showDisplayMean)
  }
  if (debugFinder)
  {
    print("here")
  }

  if (saveFile)
  {
    # plot uniform CI for quantile
    for (i in 1:unlist(optset["dDim"]))
    {
      filename <- paste(paste("QTE",toString(i)),".eps")
      if (optset["discY"] == FALSE)
      {
        plotQuantile(lenQ=length(Q)/unlist(optset["dDim"]),Qmat=Qmat[(1+(length(Q)/unlist(optset["dDim"]))*(i-1)):((length(Q)/unlist(optset["dDim"]))*i),],optset=optset,alpha=alpha,n=sum((df$T==1)&(df$X==i)),filename=filename)
      } else {
        plotQuantile(lenQ=length(QUB)/unlist(optset["dDim"]),Qmat=QUBmat[(1+(length(QUB)/unlist(optset["dDim"]))*(i-1)):((length(QUB)/unlist(optset["dDim"]))*i),],optset=optset,alpha=alpha,n=sum((df$T==1)&(df$X==i)),filename=filename)
        plotQuantile(lenQ=length(QLB)/unlist(optset["dDim"]),Qmat=QLBmat[(1+(length(QLB)/unlist(optset["dDim"]))*(i-1)):((length(QLB)/unlist(optset["dDim"]))*i),],optset=optset,alpha=alpha,n=sum((df$T==1)&(df$X==i)),filename=filename)
      }

      # filename <- paste(paste("QTE",toString(i)),"90.eps")
      # if (optset["discY"] == FALSE)
      # {
      #   plotQuantile(lenQ=length(Q)/unlist(optset["dDim"]),Qmat=Qmat[(1+(length(Q)/unlist(optset["dDim"]))*(i-1)):((length(Q)/unlist(optset["dDim"]))*i),],optset=optset,alpha=alpha*2,n=sum((df$T==1)&(df$X==i)),filename=filename)
      # } else {
      #   plotQuantile(lenQ=length(QUB)/unlist(optset["dDim"]),Qmat=QUBmat[(1+(length(QUB)/unlist(optset["dDim"]))*(i-1)):((length(QUB)/unlist(optset["dDim"]))*i),],optset=optset,alpha=alpha*2,n=sum((df$T==1)&(df$X==i)),filename=filename)
      #   plotQuantile(lenQ=length(QLB)/unlist(optset["dDim"]),Qmat=QLBmat[(1+(length(QLB)/unlist(optset["dDim"]))*(i-1)):((length(QLB)/unlist(optset["dDim"]))*i),],optset=optset,alpha=alpha*2,n=sum((df$T==1)&(df$X==i)),filename=filename)
      # }
    }
  }
  if (optset["discY"] == FALSE)
  {
    #First entry is length, and followed by CI for means
    vectorCIs <- meanCIs
    #Entry to add is the next to length of means + 1, i.e., meanCIs[1]+2
    ## length is given by: length(Q/2)*2*dDim = length(Q)*dDim
    vectorCIs[meanCIs[1]+2] <- length(Qmat[,1])*dDim
    indexStart <- meanCIs[1]+3
    # do for all x valuse
    for (xtp in 1:dDim)
    {
      # get the CI
      CI <- qCI(lenQ=length(Q)/2,Qmat=Qmat[(1+(length(Q)/2*(xtp-1))):(length(Q)/2*xtp),],nB=nBoot,alpha=alpha,n=sum((df$T==1)&(df$X==xtp)))
      #print(CI)
      # upper bound entry
      vectorCIs[indexStart:(indexStart + vectorCIs[meanCIs[1]+2]/(2*dDim) - 1)] <- CI[1,]
      # update the index
      indexStart <- indexStart + vectorCIs[meanCIs[1]+2]/(2*dDim)
      # lower bound entry
      vectorCIs[(indexStart):(indexStart + vectorCIs[meanCIs[1]+2]/(2*dDim) - 1)] <- CI[3,]
      indexStart <- indexStart + vectorCIs[meanCIs[1]+2]/(2*dDim)
      #print(vectorCIs)
    }
    return(vectorCIs)
  } else {
    return(QUBmat)
  }
}

#' (internal) Transoformation of the data.frame according to the formula specification
#'
#' This is an internal function extracting the variable names from the formula to modify the data.frame given
#' so that the formatted variable names are used in the code inside.
#'
#' @param fml a formula object
#' @param dataOrigin a data.frame object, original dataset
#' @param treatment a string, containing treatment variable varname string
#' @param posttreat a string, containing post-treatment variable varname string
#' @param cluster a string, containing cluster variable varname. Default is "" which means either its mean process or not clustered.
#'
#' @return a data.frame which contains covariates W, treatment T, posttreat D, and the outcome Y as well as the cluster.
reformatFormula <- function(fml,dataOrigin,treatment,posttreat,cluster="")
{
  fml <- Formula::as.Formula(fml)
  contVarLabels <- attr(terms(fml,rhs=1),"term.labels")
  discVarLabels <- attr(terms(fml,rhs=2),"term.labels")

  charBind <- NULL
  if (length(contVarLabels))
  {
    for (i in 1:length(contVarLabels))
    {
      varName <- paste("WCont",toString(i),sep="")
      assign(varName,pnorm((dataOrigin[,contVarLabels[i]] - mean(dataOrigin[,contVarLabels[i]]))/sqrt(var(dataOrigin[,contVarLabels[i]]))))
      if (is.null(charBind))
      {
        charBind <- paste("cbind(",varName,sep="")
      } else {
        charBind <- paste(charBind,varName,sep=",")
      }
    }
  }
  if (length(discVarLabels))
  {
    for (i in 1:length(discVarLabels))
    {
      varName <- paste("WDisc",toString(i),sep="")
      assign(varName,factor(dataOrigin[,discVarLabels[i]]))
      if (is.null(charBind))
      {
        charBind <- paste("cbind(",varName,sep="")
      } else {
        charBind <- paste(charBind,varName,sep=",")
      }
    }
  }
  charBind <- paste(charBind,")",sep="")
  W <- eval(parse(text=charBind))
  dataAdd <- data.frame(W = W)
  data <- cbind(dataOrigin,dataAdd)

  #dependent variable name, extracted as all.vars(fml)[1]]
  if (cluster == "")
  {
    dataY <- data.frame(Y = data[,all.vars(fml)[1]],T = data[,treatment],X = data[,posttreat])
  } else {
    dataY <- data.frame(Y = data[,all.vars(fml)[1]],T = data[,treatment],X = data[,posttreat],cluster= data[,cluster])
  }

  data <- cbind(data,dataY)

  return(data)
}

#' (internal) returns formula according to the spline specification
#'
#' return formula expressions for spline basis used in distribution regressions
#'
#' @param opt a list containing options
#'
#' @return a formula, fml
#'
retFormula <- function(opt)
{
  numCont <- unlist(opt["numCont"])
  numDisc <- unlist(opt["numDisc"])
  numKnots <- unlist(opt["numKnots"])
  numOrder <- unlist(opt["numOrder"])

  fml <- "dumY ~ ("
  if (numDisc > 0)
  {
    for (i in 1:numDisc)
    {
      fml <- paste(paste(fml,"W.WDisc",sep = ""),toString(i),sep = "")
      if (i < numDisc)
      {
        fml <- paste(fml,"+ ",sep = " ")
      }
    }
    if (numCont > 0)
    {
      if (opt["dumAllAdditive"] == "additive")
      {
        fml <- paste(fml,") + ",sep="")
      } else {
        fml <- paste(fml,")*(",sep="")
      }
    } else {
     fml <- paste(fml,")",sep="")
    }
  }
  if (numCont > 0)
  {
    for (i in 1:numCont)
    {
      fml <- paste(paste(fml,"bs(W.WCont",sep = ""),toString(i),sep = "")
      fml <- paste(paste(fml,",Boundary.knots=c(0,1),degree=",sep = ""),numOrder,sep = "")
      if (numKnots == 1)
      {
        fml <- paste(fml,",knots=0.5)",sep = "")
      } else if (numKnots == 2)
      {
        fml <- paste(fml,",knots=c(0.33,0.66))",sep = "")
      } else if (numKnots == 3)
      {
        fml <- paste(fml,",knots=c(0.25,0.5,0.75))",sep = "")
      }
      if (i < numCont)
      {
        if (opt["additiveSpline"] == TRUE)
        {
          fml <- paste(fml," + ",sep = "")
        } else {
          fml <- paste(fml,":",sep = "")
        }
      }
    }
    if (opt["dumAllAdditive"] == "")
    {
      fml <- paste(fml,")",sep="")
    }
  }
  return(fml)
}

#' (internal) generating distribution regressions
#'
#' distribution regression routine
#'
#' @param df a data.frame
#' @param sg a string, expression for the subgroup
#' @param itr a vector, equispaced grids
#' @param opt a list containing the options
#' @param fml a formula object
#' @param link a string, link function
#'
#' @return FCFs, a list of glm results for every grids for Y to be evaluated.
#'
genDR <- function(df,sg,itr,opt,fml,link)
{
  # list of counterfactual distribution for aech knot
  FCFs <- list()

  # formula
  if (opt["typeSieve"] != "formula")
  {
    fml <- retFormula(opt = opt)
  }

  #  objects to be used for the formula
  df$subg <- eval(parse(text=sg))

  # create subgroup dataframe
  dfSub <- df[df$subg==1,]
  dfSub <- na.omit(dfSub)

  # DR for each value of grid yi
  i <- 0
  FCFs <- list()
  if (opt["optMute"] == FALSE)
  {
    str(dfSub)
  }

  for (yi in itr)
  {
    i <- i+1
    dfSub$dumY <- (dfSub$Y <= yi)
    FY <- glm(formula=fml,data=dfSub,family=binomial(link = link), control = list(maxit = 100))
    FCFs[[i]] <- FY
  }
  # matrix of coefficients at each yi
  return(FCFs)
}

#' (internal) Two-way plot of two cdfs
#'
#' Function to return quantile for a given ranking
#' Returns list of functions of chebyshev interpolated distribution coefficients
#'
#' @param itr a vector of grids for x axis
#' @param F1 a vector of values for first objects for y axis
#' @param F2 a vector of values for second objects for y axis
#' @param label1 a string, label for F1
#' @param label2 a string, label for F2
#' @param title a string, title for the plot
#'
#' @return ggplot object readily printable.
#'
genPlot2way <- function(itr,F1,F2,label1,label2,title)
{
  # Plot averaged unconditional cdf prediction: smoothed
  plotData <- data.frame(x=itr,y=F1,case=label1)
  # against empirical cdf of unconditional cdf at finer grids
  plotDataY <- rbind(plotData,data.frame(x=itr,y=F2,case=label2))
  return(ggplot(plotDataY, aes(x=x, y=y)) + geom_line(aes(color = factor(case))) + ggtitle(title))
}


#' (internal) Checking for the fit of the predictions
#'
#' proc to check fit by ploting as well as returning FY predictions
#' It internally returns the smoothed predictions as well as plots. Plots are muted in default
#'
#' @param df a data.frame
#' @param sg a string for the subgroup expression
#' @param FYCFs a list of glm results
#' @param itrY a vector of Y grids
#' @param opt a list of options
#' @param title a string, title
#'
#' @return FYPred which is the prediction of FY for each observation.
#'
chkFit <- function(df,sg,FYCFs,itrY,opt,title)
{
  # evaluate string expression of subgroup to pick
  df$subg <- eval(parse(text=sg))
  # container
  FYWpred <- matrix(1:length(itrY)*sum(df$subg==1),sum(df$subg==1),length(itrY))
  FYWMPred <- c(1:length(itrY))
  count <- 0

  # create subgroup dataframe
  dfSub <- df[df$subg==1,]
  dfSub <- na.omit(dfSub)

  for (i in itrY)
  {
    count <- count + 1
    # prediction conditional on W: n*itrY matrix
    FYWpred[,count] <- predict(FYCFs[[count]],newdata=dfSub,type="response")
    # averaged out for n
    FYWMPred[count] <- weighted.mean(FYWpred[,count])
  }

  # evaluate on finer grids
  itrI <- seq(itrY[1],itrY[length(itrY)],length=unlist(opt["nGridsFine"]))

  FYPred <- c(itrI)
  FYE <- c(itrI)
  count <- 0

  #returing interpolated predictions
  FYPred <- approx(x=itrY,y=FYWMPred,xout=itrI,method="linear",ties=ordered)$y

  # ecdf
  for (i in 1:length(itrI))
  {
    FYE[i] <- weighted.mean((dfSub$Y <= itrI[i]))
  }

  # Rearrangement
  # FYPred <- sort(FYPred)
  # plot of smoothed averaged out cdf against ecdf
  if (opt["optMute"] == FALSE)
  {
    print(genPlot2way(itr=itrI,F1=FYPred,F2=FYE,label1="smoothed",label2="ecdf",title=title))
  }
  return(FYPred)
}

#' (internal) Showing the treatment effects
#'
#' Return Treatment Effects
#'
#' Muted in default.
#'
#' @param df data.frame
#' @param Pred1 prediction of Y1
#' @param Pred0 prediction of Y0
#' @param seq1 a string, subgroup expression for Y1
#' @param seq0 a string, subgroup expression for Y0
#' @param itr grids for Y
#' @param opt a list of options
#'
#' @return None.
#'
showDiff <- function(df,Pred1,Pred0,seq1,seq0,itr,opt)
{
  itr <- seq(itr[1],itr[length(itr)],length=unlist(opt["nGridsFine"]))

  plotDataY0 <- data.frame(x=itr,y=Pred0,case=0)
  plotDataY1 <- data.frame(x=itr,y=Pred1,case=1)
  plotDataTE <- rbind(plotDataY1,plotDataY0)
  df$subg1 <- eval(parse(text=seq1))
  df$subg0 <- eval(parse(text=seq0))

  if (opt["optMute"] == FALSE)
  {
    print(genPlot2way(itr=itr,F1=Pred0,F2=Pred1,label1="Y0",label2="Y1",title="Treatment Effects"))

    print(ggplot(plotDataTE, aes(x=x, y=y)) + geom_line(aes(color = factor(case))))
    cat(mean(df$Y[df$subg1 == 1]),"\n")
    cat(mean(df$Y[df$subg0 == 1]),"\n")
    cat(mean(df$Y[df$subg1 == 1]) - mean(df$Y[df$subg0 == 1]),"\n")
  }
}

#' (internal) Prediction for the conditional cdfs
#'
#' prediction of FY conditional on W
#'
#' @param df data.frame
#' @param sg a string, expression for the subgroup
#' @param FYCFs a glm results, distribution regression coefficients
#' @param itr a vector, grids
#' @param opt a list, options
#' @param n an integer, subsample size
#'
#' @return FYWpred, prediction of FY for each observation
#'
predFY <- function(df,sg,FYCFs,itr,opt,n)
{

  FYWpred <- matrix(1:length(itr)*n,n,length(itr))
  df$subg <- eval(parse(text=sg))
  # create subgroup dataframe
  dfSub <- df[df$subg==1,]
  dfSub <- na.omit(dfSub)
  ## count levels
  # print(sapply(fctr, nlevels))
  tt <- tryCatch(predict(FYCFs[[1]],newdata=dfSub,type="response"),error=function(e) e, warning=function(w) w)
  if(is(tt,"warning"))
  {
    print(tt)
    print(sg)
    FYWpred <- "RankDef"
  } else {
    for (count in 1:length(itr))
    {
      FYWpred[,count] <- predict(FYCFs[[count]],newdata=dfSub,type="response")
    }
  }
  return(FYWpred)
}

#' (internal) Numerical integration with respect to the estimated cdf
#'
#' Numerical integration for usual grids
#'
#' @param itr grids
#' @param FY values
#' @param opt a list of options
#'
#' @return mean of integrated value with respect to the values FY
#'
integrate1 <- function(itr,FY,opt)
{
  mean <- 0
  iStart <- 1
  iEnd <- (length(itr)-1)
  #if (opt["discY"] == TRUE)
  #{
  #  mean <- mean + itr[1]*FY[1]
  #  iEnd <- (length(itr)-2)
  #}
  for (i in iStart:iEnd)
  {
    mean <- mean + itr[i]*(FY[i+1] - FY[i])
  }
  #if (opt["discY"] == FALSE)
  #{
  mean <- mean + itr[length(itr)]*(1 - FY[length(itr)])
  #}
  return(mean)
}

#' (internal) Take numerical inverse of the cdf
#'
#' get quantile by inversion
#'
#' @param FY a vector of values
#' @param kI a vector of grids for FY
#' @param kU a vector of grids to evaluate
#'
#' @return QY quantile of FY
#'
getQuantile <- function(FY,kI,kU)
{
  kUext <- c(0.01,kU)
  QY <- c(1:length(kU))*0
  FY <- sort(FY)
  curj = 1
  for (i in 1:(length(kU)+1))
  {
    for (j in curj:length(kI))
    {
      if (FY[j] >= kUext[i])
      {
        if (i > 1)
        {
          QY[i-1] <- kI[j]
        }
        curj <- j
        break
      }
    }
  }
  return(QY)
}

#' (internal) Main caller for the estimation procedure
#'
#' the common procedure called by both mainestimate and bootstrapProc
#'
#' @param main data.frame of the reformatted dataset
#' @param optset a list of options
#' @param fml a formula reformatted ready for the glm
#' @param Ytype a string, either "Y1" or "Yb"
#' @param link a string, link function for the glm
#'
#' @return reporting of the list of mean estimates, $reportingV and $reportingQ

mainProcedure <- function(main,optset,fml,Ytype,link)
{
  if (Ytype != "Y1")
  {
    # replace T = 0 of alternative Y measure with original Y
    main$Yb[main$T==0] <- main$Y[main$T == 0]
    # maintain the original Y
    Yorig <- main$Y
    # keep the original data.frame as mainOrg: used for estimating factual quantiles
    mainOrg <- main
    # this is used for computing mean
    main <- cbind(main,Yorig)
    # use Yb as Y
    main$Y <- main$Yb
    # for factual Y, we must refer to Yorig
  }

  ## getting grids
  # genearate knots
  if (optset["discY"] == FALSE)
  {
    itrY <- seq(from = min(main$Y), to = max(main$Y), length=unlist(optset["nGrids"]))
    itrI <- seq(from = min(main$Y), to = max(main$Y), length=unlist(optset["nGridsFine"]))
  } else {
    #itrY <- c(0,gridsGen(min = min(main$Y[main$Y > 0]), max = max(main$Y), nG=unlist(optset["nGrids"])))
    #itrI <- c(0,gridsGen(min = min(main$Y[main$Y > 0]), max = max(main$Y), nG=unlist(optset["nGridsFine"])))
    itrY <- seq(from = min(main$Y), to = max(main$Y), length=unlist(optset["nGrids"]))
    itrI <- seq(from = min(main$Y), to = max(main$Y), length=unlist(optset["nGridsFine"]))
  }
  itrU <- seq(from = 0, to = 1, length=unlist(optset["nGridsFine"]))

  # Common parts
  n0 <- sum(main$T==0)
  n1 <- sum(main$T==1)
  n <- n0 + n1
  ## getting the distribution regressions
  # ith column is ith grids coefficients
  # Dist Regression for Y0
  optset0 <- optset
  #optset0["numKnots"] <- 2
  #optset0["numOrder"] <- 2
  FY0CFs <-genDR(df=main,sg="(df$T==0)",itr=itrY,opt=optset0,fml=fml,link=link)
  # Dist Regression for Y1
  FY1CFs <-genDR(df=main,sg="(df$T==1)",itr=itrY,opt=optset0,fml=fml,link=link)
  ## predictions based on estimated DRs
  FY0Pred <- chkFit(df = main,sg = "(df$T==0)",FYCFs = FY0CFs,itrY = itrY,opt = optset0,title="Y0")
  FY1Pred <- chkFit(df = main,sg = "(df$T==1)",FYCFs = FY1CFs,itrY = itrY,opt = optset0,title="Y1")
  ## ATE estimated
  if (optset["optMute"] == FALSE)
  {
    showDiff(df=main,Pred0=FY0Pred,Pred1=FY1Pred,seq1="(df$T==1)",seq0="(df$T==0)",itr=itrY,opt=optset)
  }

  if (Ytype == "Y1")
  {
    mean1 <- weighted.mean(main$Y[main$T==1])
    mean0 <- weighted.mean(main$Y[main$T==0])
  } else {
    mean1 <- weighted.mean(main$Yorig[main$T==1])
    mean0 <- weighted.mean(main$Yorig[main$T==0])
  }
  itrUShow <- seq(unlist(optset["plotBeg"]),unlist(optset["plotEnd"]),by=unlist(optset["plotBy"]))
  reportFlag <- 0
  for (i in 1:unlist(optset["dDim"]))
  {
    # sample size variables updated
    nX <- sum((main$T==1)&(main$X==i))
    sgExp <- paste(paste("(df$T==1)&(df$X == ",toString(i)),")")
    # Dist Regression for Y1X1
    optsetX <- optset
    FY1XCFs <-genDR(df=main,sg=sgExp,itr=itrY,opt=optset,fml=fml,link=link)
    FY1XPred <- chkFit(df = main,sg = sgExp,FYCFs = FY1XCFs,itrY = itrY,opt = optset,title="Y1X")
    # FY predicted
    FY0Wpred1 <- predFY(df=main,sg=sgExp,FYCFs=FY0CFs,itr=itrY,opt=optset,n=nX)
    if (FY0Wpred1 == "RankDef")
    {
       reportFlag <- NA
    }
    FY1Wpred1 <- predFY(df=main,sg=sgExp,FYCFs=FY1CFs,itr=itrY,opt=optset,n=nX)
    if (FY1Wpred1 == "RankDef")
    {
      reportFlag <- NA
    }
    FY1XWpred1 <- predFY(df=main,sg=sgExp,FYCFs=FY1XCFs,itr=itrY,opt=optset,n=nX)
    if (FY1XWpred1 == "RankDef")
    {
      reportFlag <- NA
    }
    if (!is.na(reportFlag))
    {
      ## Predicting quantiles
      # Cpp alternative
      if (optset["discY"] == FALSE)
      {
        QY1XCF <- predQ1Cpp(FY1Wpred1,FY0Wpred1,itrY,itrU,unlist(optset["nGridsFine"]),nX,0,0)
      } else {
        QY1XCFUB <- predQ1Cpp(FY1Wpred1,FY0Wpred1,itrY,itrU,unlist(optset["nGridsFine"]),nX,1,1)
        QY1XCFLB <- predQ1Cpp(FY1Wpred1,FY0Wpred1,itrY,itrU,unlist(optset["nGridsFine"]),nX,1,0)
      }
      ## Evaluate counterfactual distributions
      if (optset["discY"] == FALSE)
      {
        FY0XCFM <- predCFCpp(FY1XWpred1,QY1XCF,itrY,itrU,c(1:nX)*0,unlist(optset["nGridsFine"]),nX,0,0)
      } else {
        FY0XCFMUB <- predCFCpp(FY1XWpred1,QY1XCFUB,itrY,itrU,c(1:nX)*0,unlist(optset["nGridsFine"]),nX,0,1)
        FY0XCFMLB <- predCFCpp(FY1XWpred1,QY1XCFLB,itrY,itrU,c(1:nX)*0,unlist(optset["nGridsFine"]),nX,0,1)
      }

      if (Ytype != "Y1")
      {
        FY1XCFsOrg <-genDR(df=mainOrg,sg=sgExp,itr=itrY,opt=optset,fml=fml,link=link)
        FY1XPredOrg <- chkFit(df = mainOrg,sg = sgExp,FYCFs = FY1XCFsOrg,itrY = itrY,opt = optset,title="Y1X")
      }

      if (optset["discY"] == FALSE)
      {
        QY0XCFM <- getQuantileCpp(FY0XCFM,itrI,itrUShow)
        if (Ytype == "Y1")
        {
          QY1XPred <- getQuantileCpp(FY1XPred,itrI,itrUShow)
        } else {
          QY1XPred <- getQuantileCpp(FY1XPredOrg,itrI,itrUShow)
        }
      } else {
        QY0XCFMUB <- getQuantileCpp(FY0XCFMUB,itrI,itrUShow)
        QY1XPred <- getQuantileCpp(FY1XPred,itrI,itrUShow)
        QY0XCFMLB <- getQuantileCpp(FY0XCFMLB,itrI,itrUShow)
      }

      if (i == 1)
      {
        if (optset["discY"] == FALSE)
        {
          meanC0 <- integrate1(itrI,FY0XCFM,opt=optset)
        } else {
          meanC0UB <- integrate1(itrI,FY0XCFMUB,opt=optset)
          meanC0LB <- integrate1(itrI,FY0XCFMLB,opt=optset)
        }
        if (Ytype == "Y1")
        {
          meanE1 <- weighted.mean(main$Y[(main$T==1)&(main$X==i)])
        } else {
          meanE1 <- weighted.mean(main$Yorig[(main$T==1)&(main$X==i)])
        }
        if (optset["discY"] == FALSE)
        {
          reportingQ <- (QY1XPred - QY0XCFM)
          reportingV <- c(meanE1 - meanC0)
          ATEPred <- (meanE1 - meanC0)*nX/n1
        } else {
          reportingQUB <- (QY1XPred - QY0XCFMUB)
          reportingVUB <- c(meanE1 - meanC0UB)
          ATEPredUB <- (meanE1 - meanC0UB)*nX/n1
          reportingQLB <- (QY1XPred - QY0XCFMLB)
          reportingVLB <- c(meanE1 - meanC0LB)
          ATEPredLB <- (meanE1 - meanC0LB)*nX/n1
        }
      } else
      {
        if (optset["discY"] == FALSE)
        {
          meanC0 <- integrate1(itrI,FY0XCFM,opt=optset)
        } else {
          meanC0UB <- integrate1(itrI,FY0XCFMUB,opt=optset)
          meanC0LB <- integrate1(itrI,FY0XCFMLB,opt=optset)
        }

        if (Ytype == "Y1")
        {
          meanE1 <- weighted.mean(main$Y[(main$T==1)&(main$X==i)])
        } else {
          meanE1 <- weighted.mean(main$Yorig[(main$T==1)&(main$X==i)])
        }

        if (optset["discY"] == FALSE)
        {
          reportingQ <- c(reportingQ,QY1XPred - QY0XCFM)
          reportingV <- c(reportingV,meanE1 - meanC0)
          ATEPred <- ATEPred + (meanE1 - meanC0)*nX/n1
        } else {
          reportingQUB <- c(reportingQUB,QY1XPred - QY0XCFMUB)
          reportingVUB <- c(reportingVUB,meanE1 - meanC0UB)
          ATEPredUB <- ATEPredUB + (meanE1 - meanC0UB)*nX/n1
          reportingQLB <- c(reportingQLB,QY1XPred - QY0XCFMLB)
          reportingVLB <- c(reportingVLB,meanE1 - meanC0LB)
          ATEPredLB <- ATEPredLB + (meanE1 - meanC0LB)*nX/n1
        }
      }
    } else {
        reportingV <- seq(length = unlist(optset["dDim"]))*NA
        reportingQ <- seq(length = itrUShow*unlist(optset["dDim"]))*NA
        print("Rank Deficiency detected for this dataset. This estimation will not be used.")
    }
  }

  if (optset["discY"] == FALSE)
  {
    reportingV <- c(reportingV,mean1 - mean0,ATEPred)
    reporting <- list("reportingV"=reportingV,"reportingQ"=reportingQ)
  } else {
    reportingVUB <- c(reportingVUB,mean1 - mean0,ATEPredUB)
    reportingVLB <- c(reportingVLB,mean1 - mean0,ATEPredLB)
    reporting <- list("reportingVUB"=reportingVUB,"reportingVLB"=reportingVLB,"reportingQUB"=reportingQUB,"reportingQLB"=reportingQLB)
  }
  return(reporting)
}


#' (internal) Reformatting bootstrapped quantile estimates into the CIs
#'
#' a function reformat matrix of quantile estimates into the uniform CIs
#'
#' @param lenQ an integer, length of Q
#' @param Qmat a matrix, matrix of bootstrapped quantile estimates
#' @param alpha a float, the size of the test
#' @param nB an integer, number of bootstrap iterations
#' @param n an integer, sample size
#'
#' @return a matrix of three columns of CIs with mean estimates in the middle.
#'
qCI <- function(lenQ,Qmat,alpha,nB,n)
{
  iqrange <- c(1:lenQ)
  tstat <- matrix(lenQ*(nB-1),lenQ,nB-1)
  for (i in 1:lenQ)
  {
    Z <- sqrt(n)*(Qmat[i,2:nB]-Qmat[i,1])
    iqrange[i] <- (quantile(Z,0.75,na.rm=TRUE) - quantile(Z,0.25,na.rm=TRUE))/(qnorm(0.75) - qnorm(0.25))
    tstat[i,] <- sqrt(n)*abs(Qmat[i,2:nB]-Qmat[i,1])/(iqrange[i])
  }
  tstatH <- c(1:(nB-1))
  for (b in 1:(nB-1))
  {
    tstatH[b] <- max(tstat[,b])
  }
  CI <- matrix(lenQ*3,3,lenQ)

  for (i in 1:lenQ)
  {
    CI[1,i] <- Qmat[i,1] - quantile(tstatH,1-alpha,na.rm=TRUE)*iqrange[i]/sqrt(n)
    CI[2,i] <- Qmat[i,1]
    CI[3,i] <- Qmat[i,1] + quantile(tstatH,1-alpha,na.rm=TRUE)*iqrange[i]/sqrt(n)
  }
  return(CI)
}

#' (internal) Generating and saving the plots of the quantiles with uniform CIs
#' generating plot of quantiles with uniform CI
#'
#' @param lenQ an integer, length of Q
#' @param Qmat a matrix, matrix of bootstrapped quantile estimates
#' @param optset a list of options
#' @param alpha a float, the size of the test
#' @param n an integer, sample size
#' @param filename a string for the path to save the plot figures
#'
#' @return none.

plotQuantile <- function(lenQ,Qmat,optset,alpha,n,filename)
{
  CI <- qCI(lenQ=lenQ,Qmat=Qmat,nB=unlist(optset["nBoot"]),alpha=alpha,n=n)
  plotDt <- data.frame(Quantile = seq(unlist(optset["plotBeg"]),unlist(optset["plotEnd"]),by=unlist(optset["plotBy"])),TE=CI[1,],plots="LB")
  plotDt <- rbind(plotDt,data.frame(Quantile = seq(unlist(optset["plotBeg"]),unlist(optset["plotEnd"]),by=unlist(optset["plotBy"])),TE=CI[2,],plots="Estimates"))
  plotDt <- rbind(plotDt,data.frame(Quantile = seq(unlist(optset["plotBeg"]),unlist(optset["plotEnd"]),by=unlist(optset["plotBy"])),TE=CI[3,],plots="UB"))
  # ggplot(plotDt,aes(x=x,y=y)) + geom_line(aes(color=factor(case)))
  print(ggplot(plotDt,aes(x=Quantile,y=TE)) + geom_line(aes(linetype=plots),color="#099191") + scale_linetype_manual(values=c("dotted","solid","dotted")) + theme(legend.position="top"))
  ggsave(filename,plot = last_plot())
}

#' (internal) Displaying the mean estimation tables with CIs
#' generating table of mean estimates with CIs
#'
#' @param repV reportingV object
#' @param optset a list of options
#' @param n the sample size
#' @param repVLB reportingV object for the partial identification. Default is NA
#' @param curB an integer, current value of the bootstrap procedure
#' @param showDisplayMean a boolean, if TRUE, then it shows an intermediary table
#'
#' @return None.
#'
displayMeans <- function(repV,optset,n,repVLB=NA,curB,showDisplayMean=TRUE)
{
  # rescaled interquartile range
  lenVector <- (unlist(optset["dDim"])+2)
  iqrange <- c(1:lenVector)
  tstat <- matrix(lenVector*(curB-1),lenVector,curB-1)
  iqrangeLB <- c(1:lenVector)
  tstatLB <- matrix(lenVector*(curB-1),lenVector,curB-1)
  for (i in 1:lenVector)
  {
    Z <- sqrt(n)*(repV[i,2:curB]-repV[i,1])
    iqrange[i] <- (quantile(Z,0.75,na.rm=TRUE) - quantile(Z,0.25,na.rm=TRUE))/(qnorm(0.75) - qnorm(0.25))
    tstat[i,] <- sqrt(n)*abs(repV[i,2:curB] - repV[i,1])/(iqrange[i])
    if (optset["discY"] == TRUE)
    {
      iqrangeLB[i] <- sqrt(n)*(quantile(repVLB[i,2:curB] - repVLB[i,1],0.75,na.rm=TRUE) - quantile(repVLB[i,2:curB] - repVLB[i,1],0.25,na.rm=TRUE))/(qnorm(0.75) - qnorm(0.25))
      tstatLB[i,] <- sqrt(n)*abs(repVLB[i,2:curB] - repVLB[i,1])/(iqrange[i])
    }
  }

  if (optset["discY"] == TRUE)
  {
    message("Upper bounds \n")
    cat("Upper bounds \n")
  }
  if (showDisplayMean)
  {
    cat(sprintf("%6s \t \t %6s \t %6s \t %6s \t %6s \t %6s \t %6s \n","","mean","se","(95%","CI)","(90%","CI)"))
    for (i in 1:lenVector)
    {
      if (i <= unlist(optset["dDim"]))
      {
        strName <- paste("CATE:",toString(i),"\t")
      }
      if (i == unlist(optset["dDim"]) + 1)
      {
        strName <- "ATE \t \t"
      }
      if (i == unlist(optset["dDim"]) + 2)
      {
        strName <- "ATEPred \t"
      }
      cat("sd based \n")
      cat(sprintf("%6s \t %6.4f \t %6.4f \t (%6.4f \t %6.4f) \t (%6.4f  \t %6.4f) \n",strName,repV[i,1],sd(repV[i,2:curB] - repV[i,1],na.rm=TRUE),
                  repV[i,1] - sd(repV[i,2:curB] - repV[i,1],na.rm=TRUE)*1.96,
                  repV[i,1] + sd(repV[i,2:curB] - repV[i,1],na.rm=TRUE)*1.96,
                  repV[i,1] - sd(repV[i,2:curB] - repV[i,1],na.rm=TRUE)*1.645,
                  repV[i,1] + sd(repV[i,2:curB] - repV[i,1],na.rm=TRUE)*1.645))
    }
  }
  # for (i in 1:lenVector)
  # {
  #   if (i <= unlist(optset["dDim"]))
  #   {
  #     strName <- paste("CATE:",toString(i),"\t")
  #   }
  #   if (i == unlist(optset["dDim"]) + 1)
  #   {
  #     strName <- "ATE \t \t"
  #   }
  #   if (i == unlist(optset["dDim"]) + 2)
  #   {
  #     strName <- "ATEPred \t"
  #   }
  #   message("quantile based \n")
  #   message(sprintf("%6s \t %6.4f \t %6.4f \t (%6.4f \t %6.4f) \t (%6.4f  \t %6.4f) \n",strName,repV[i,1],sd(repV[i,2:curB] - repV[i,1],na.rm=TRUE),
  #                   repV[i,1] - quantile(tstat[i,],1-0.05,na.rm=TRUE)*iqrange[i]/sqrt(n),
  #                   repV[i,1] + quantile(tstat[i,],1-0.05,na.rm=TRUE)*iqrange[i]/sqrt(n),
  #                   repV[i,1] - quantile(tstat[i,],1-0.1,na.rm=TRUE)*iqrange[i]/sqrt(n),
  #                   repV[i,1] + quantile(tstat[i,],1-0.1,na.rm=TRUE)*iqrange[i]/sqrt(n)))
  #
  #   cat("quantile based \n")
  #   cat(sprintf("%6s \t %6.4f \t %6.4f \t (%6.4f \t %6.4f) \t (%6.4f  \t %6.4f) \n",strName,repV[i,1],sd(repV[i,2:curB] - repV[i,1],na.rm=TRUE),
  #               repV[i,1] - quantile(tstat[i,],1-0.05,na.rm=TRUE)*iqrange[i]/sqrt(n),
  #               repV[i,1] + quantile(tstat[i,],1-0.05,na.rm=TRUE)*iqrange[i]/sqrt(n),
  #               repV[i,1] - quantile(tstat[i,],1-0.1,na.rm=TRUE)*iqrange[i]/sqrt(n),
  #               repV[i,1] + quantile(tstat[i,],1-0.1,na.rm=TRUE)*iqrange[i]/sqrt(n)))
  # }

  #Var-Cov matrix of CATE1 and CATE2
  # S <- var(cbind(repV[1,2:curB],repV[2,2:curB]),na.rm=TRUE)
  # if ((curB == unlist(optset["nBoot"])) & (unlist(optset["nBoot"]) >= 100))
  # {
  #   numGrids <- 100
  #   s2Range <- seq(-1,2,length=numGrids)
  #   s1Range <- seq(-0.5,1.5,length=numGrids)
  #   CR <- data.frame()
  #   for (s1 in 1:numGrids)
  #   {
  #     for (s2 in 1:numGrids)
  #     {
  #       Fstat <- (cbind(repV[1,1] - s1Range[s1],repV[2,1] - s2Range[s2]))%*%solve(S)%*%t(cbind(repV[1,1] - s1Range[s1],repV[2,1] - s2Range[s2]))
  #       X1 <- s1Range[s1]
  #       X2 <- s2Range[s2]
  #       if (pchisq(Fstat,df=2) <= 0.95)
  #       {
  #         CRnewRow <- data.frame(X1 = X1, X2 = X2,r=1)
  #       } else         {
  #         CRnewRow <- data.frame(X1 = X1, X2 = X2,r=0)
  #       }
  #       CR <- rbind(CR,CRnewRow)
  #     }
  #   }
  #   #print(CR)
  #   print(ggplot(CR,aes(x = X1, y = X2)) + geom_point(aes(color = factor(r)),shape=18))
  # }
  # cat("Covariance in CATEs \t %6.4f \n", var(repV[1,2:curB],repV[2,2:curB],na.rm = TRUE))
  if (optset["discY"] == TRUE)
  {
    message("----------------------------------------- \n")
    message("Lower bounds \n")
    cat("----------------------------------------- \n")
    cat("Lower bounds \n")
    for (i in 1:lenVector)
    {
      if (i <= unlist(optset["dDim"]))
      {
        strName <- paste("CATE:",toString(i),"\t")
      }
      if (i == unlist(optset["dDim"]) + 1)
      {
        strName <- "ATE \t \t"
      }
      if (i == unlist(optset["dDim"]) + 2)
      {
        strName <- "ATEPred \t"
      }
      message(sprintf("%6s \t %6.4f \t %6.4f \t (%6.4f \t %6.4f) \t (%6.4f  \t %6.4f) \n",strName,repVLB[i,1],sd(repVLB[i,2:curB] - repVLB[i,1],na.rm=TRUE),
                  repVLB[i,1] - sd(repVLB[i,2:curB] - repVLB[i,1],na.rm=TRUE)*1.96,
                  repVLB[i,1] + sd(repVLB[i,2:curB] - repVLB[i,1],na.rm=TRUE)*1.96,
                  repVLB[i,1] - sd(repVLB[i,2:curB] - repVLB[i,1],na.rm=TRUE)*1.645,
                  repVLB[i,1] + sd(repVLB[i,2:curB] - repVLB[i,1],na.rm=TRUE)*1.645))

      cat(sprintf("%6s \t %6.4f \t %6.4f \t (%6.4f \t %6.4f) \t (%6.4f  \t %6.4f) \n",strName,repVLB[i,1],sd(repVLB[i,2:curB] - repVLB[i,1],na.rm=TRUE),
                  repVLB[i,1] - sd(repVLB[i,2:curB] - repVLB[i,1],na.rm=TRUE)*1.96,
                  repVLB[i,1] + sd(repVLB[i,2:curB] - repVLB[i,1],na.rm=TRUE)*1.96,
                  repVLB[i,1] - sd(repVLB[i,2:curB] - repVLB[i,1],na.rm=TRUE)*1.645,
                  repVLB[i,1] + sd(repVLB[i,2:curB] - repVLB[i,1],na.rm=TRUE)*1.645))
    }
  }

  if (optset["discY"] == TRUE)
  {
    message("----------------------------------------- \n")
    cat("----------------------------------------- \n")
    message(sprintf("%6s \t \t %6s \t %6s \t %6s \t %6s \t %6s \t %6s \n","","mean LB","mean UB","(95%","CI)","(90%","CI)"))
    cat(sprintf("%6s \t \t %6s \t %6s \t %6s \t %6s \t %6s \t %6s \n","","mean LB","mean UB","(95%","CI)","(90%","CI)"))
    for (i in 1:lenVector)
    {
      if (i <= unlist(optset["dDim"]))
      {
        strName <- paste("CATE:",toString(i),"\t")
      }
      if (i == unlist(optset["dDim"]) + 1)
      {
        strName <- "ATE \t \t"
      }
      if (i == unlist(optset["dDim"]) + 2)
      {
        strName <- "ATEPred \t"
      }
      message(sprintf("%6s \t %6.4f \t %6.4f \t (%6.4f \t %6.4f) \t (%6.4f  \t %6.4f) \n",strName,repVLB[i,1],repV[i,1],
                  quantile(repVLB[i,2:curB],0.025,na.rm=TRUE),
                  quantile(repV[i,2:curB],0.975,na.rm=TRUE),
                  quantile(repVLB[i,2:curB],0.05,na.rm=TRUE),
                  quantile(repV[i,2:curB],0.95,na.rm=TRUE)))
      cat(sprintf("%6s \t %6.4f \t %6.4f \t (%6.4f \t %6.4f) \t (%6.4f  \t %6.4f) \n",strName,repVLB[i,1],repV[i,1],
                  quantile(repVLB[i,2:curB],0.025,na.rm=TRUE),
                  quantile(repV[i,2:curB],0.975,na.rm=TRUE),
                  quantile(repVLB[i,2:curB],0.05,na.rm=TRUE),
                  quantile(repV[i,2:curB],0.95,na.rm=TRUE)))
    }
  }

  meanCIs <- c()
  #first entry contains the length of means
  meanCIs <- lenVector*2
  for (i in 1:lenVector)
  {
    meanCIs[2*i] <- repV[i,1] - sd(repV[i,2:curB] - repV[i,1],na.rm=TRUE)*1.96
    meanCIs[2*i+1] <- repV[i,1] + sd(repV[i,2:curB] - repV[i,1],na.rm=TRUE)*1.96
  }
  return(meanCIs)
}

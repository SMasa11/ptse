# ptse
R package for Post-Treatment Subgroup Effect analysis

## Instal
```devtools::install_github("SMasa11/ptse")```

## Usage
The function `ptsEst()` returns both
1. Estimation of average subgroup effects and quantile subgroup differences
2. Bootstrap Confidence Intervals.

### Example
```
CIs <- ptsEst(df=main,dDim=2,treatment="treatment",posttreat="client",
clusterInference=TRUE,cluster="cluster",numKnots=1,numOrder=1,formula=formula,Ytype="Yb")
```
This procedure returns the list of CIs, while mean estimates, tables for CIs and quantile plots will be produced during the procedures.

#### dDim: the number of the subgroups
For example, dDim = 2 is the case where the subgroup is defined by a binary "posttreat" variable.

#### formula: formula object which should have the form
 outcome ~ continuousVariables | discreteVariables
 continuousVariables and discreteVariables may be multiple variables combined with "+"
   example: Y ~ cont1 + cont2 + cont3 | disc1 + disc2
 If there is no continuous variable, or discrete variable, replace it with "0"
   example: Y ~ cont1 + cont2 | 0
   example: Y ~ 0 | disc1 + disc2

#### numKnots and numOrder specify the number of knots and the order of the spline basis
 spline basis will be generated for continuousVariables only
 where the spline basis of continuous variables will be multiplied with (discreteVariables)
  Expected format will be, splineBasis*(disc1 + disc2 + disc3)
  remark: if the cross terms with discrete variables will be generated "as is".
   If you wish to allow splines to differe for all the combinations of discrete variables, specify so in the discreteVariables part of the "formula".

#### df: data.frame object which contains
 treatment: a variable specified as the argument treatment, e.g., "treatment" (T in the paper)
 posttreat: a variable specified as the argument posttreat, e.g., "client" (D in the paper)
 outcome: a variable specified as the left-hand side of the "formula", e.g. "output" (Y in the paper)
 Yb: the proxy variable, the proxy variable must be named as "Yb" in your data.frame
 cluster: a variable speficied as the argument cluster, e.g., "cluster"
  If clusterInference = FALSE, then cluster argument must be "" 

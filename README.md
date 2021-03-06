# ptse
R package for Post-Treatment Subgroup Effect analysis, Sawada (2018)

 For the analysis, refer the paper "Identification and inference of post-treatment subgroup effects." on my website: https://sites.google.com/view/masayukisawada/research

## Instal
```devtools::install_github("SMasa11/ptse")```

## Usage
The function `ptsEst()` returns both
1. Estimation of average subgroup effects and quantile subgroup differences
2. Bootstrap Confidence Intervals.

### Example
```
CIs <- ptsEst(df=main,dDim=2,treatment="treatment",posttreat="client",
clusterInference=TRUE,cluster="cluster",numKnots=2,numOrder=2,formula=formula,Ytype="Yb")
```
This procedure returns the list of CIs. 

Mean estimates, tables for CIs and quantile plots will be produced during the procedures.

#### dDim: the number of the subgroups
For example, dDim = 2 is the case where the subgroup is defined by a binary "posttreat" variable.

#### formula: formula object which should have the form
 formula must take the following format

 ```Y ~ continuousVariables | discreteVariables```
 
 where continuousVariables and discreteVariables may be multiple variables combined with "+"
 
 Example :
 ``` Y ~ cont1 + cont2 + cont3 | disc1 + disc2 ```
 
 If there is no continuous variable, or discrete variable, replace it with "0"
 ``` 
     # when there is no discrete covariates
     Y ~ cont1 + cont2 | 0
     
     # when there is no continuous covariates
     Y ~ 0 | disc1 + disc2
 ```

#### numKnots and numOrder specify the number of knots and the order of the spline basis
 Spline basis will be generated for continuousVariables only
 where the spline basis of continuous variables will be multiplied with (discreteVariables). The example above means the spline basis with two knots and up to squared terms.
 
  Expected format will be, splineBasis*(disc1 + disc2 + disc3)
  
  
  remark: if the cross terms with discrete variables will be generated "as is".
  
  If you wish to allow splines to differe for all the combinations of discrete variables, specify so in the discreteVariables part of the "formula".

#### df: data.frame object which contains
 treatment: a variable specified as the argument treatment, e.g., "treatment" (T in the paper)
 
 posttreat: a variable specified as the argument posttreat, e.g., "client" (D in the paper)
 
 outcome: a variable specified as the left-hand side of the "formula", e.g. "Y" (Y in the paper)
 
 Yb: the proxy variable, the proxy variable must be named as "Yb" in your data.frame
 
 cluster: a variable speficied as the argument cluster, e.g., "cluster".
  If clusterInference = FALSE, then cluster argument must be "" 

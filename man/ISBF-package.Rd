\name{ISBF-package}
\alias{ISBF-package}
\alias{ISBF}
\docType{package}
\title{
Iterative Selection of Blocks of Features - ISBF
}
\description{
Selection of features for sparse regression estimation (like the LASSO). Selection of blocks of features when the regression parameter is sparse and constant by blocks (like the Fused-LASSO). Application to CGH data. TEST VERSION.
}
\details{
\tabular{ll}{
Package: \tab ISBF\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-11-10\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
The package contains the following functions:

1) isbfReg - performs regression estimation in the model Y = Xb + e where b is sparse.
> isbfReg(X,Y)
If b is also constant by block, the function may find blocks up to a given size K.
> isbf(X,Y,K=...)

2) isbf - particular case where X is the identity matrix, Y = b + e, and b is sparse. This function is much faster than isbfReg.
> isbf(Y)
If b is also constant by block, the function may find blocks up to size K.
> isbf(Y,K=...)

[Be careful, for functions 1 and 2, the computation time and the memory used grows with K!!!!]

3) cghISBF - applies isfb to every chromosome in a cgh array. The object returned has the type cghFLasso used in the package cghFLasso
(see Tibshirani and Wang, 2008). Therefore, it can be plotted using the package cghFLasso.

Finally, CGHDisease1 is an example of cgh array taken from the cghFLasso package.
}
\author{
Pierre Alquier <alquier@ensae.fr>
}
\references{
P. Alquier, An Algorithm for Iterative Selection of Blocks of Features, Proceedings of ALT'10, 2010, M. Hutter, F. Stephan, V. Vovk and T. Zeugmann Eds., Lecture Notes in Artificial Intelligence, pp. 35-49, Springer.

P. Alquier, Iterative Feature Selection in Least Square Regression Estimation, Annales de l'IHP, B (Proba. Stat.), 2008, vol. 44, no. 1, pp 47-88.

R. Tibshirani and P. Wang, Spatial Smoothing and Hot Spot Detection for CGH Data Using the Fused Lasso, Biostatistics, 2008, vol. 9, no. 1, pp 18-29.
}

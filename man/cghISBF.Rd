\name{cghISBF}
\alias{cghISBF}
\title{Iterative Selection of Blocks of Features for CGH arrays}
\description{
A function to process reference CGH arrays using the isbf function. Note that this function
takes the same input than the function cghFLasso of the package cghFLasso, see Tibshirani
and Wang (2008). The output is an object of the type cghFLasso. Therefore, it can be plotted
using the package cghFLasso.
}
\usage{
cghISBF(CGH.Array, chromosome, nucleotide.position, epsilon = 0.05,
K = 1, impmin = 1/100, s = NULL, v = NULL)
}
\arguments{
  \item{CGH.Array}{Numeric vector. The result of one or multiple CGH experiments.
 Each column is the log2 ratios returned from one array experiment and is ordered
according to the gene/clones position on the genome. No missing values allowed.}
  \item{chromosome}{Numeric vector. Length should be the same as CGH.Array.
The chromosome number of each gene/clone.}
  \item{nucleotide.position}{Numeric vector. Length should be the same as CGH.Array.
The nucleotide position of each gene/clone. This information is used mainly for plot.
}
  \item{epsilon}{The confidence level when running ISBF. The theoretical guarantees in Alquier (2010) is that each iteration of the ISBF procedure gets closer to the real parameter b with probability at least 1-epsilon. When epsilon is very small, the procedure becomes very conservative. When epsilon is too large, there is a risk of overfitting. If not specified, epsilon = 5\%.
}
  \item{K}{The maximal length of blocks checked in the iterations. If not specified, K=1. If K is larger than the length of a chromosome, then it
will we adjusted when the function isbf is called on this chromosome.
}
  \item{impmin}{Criterion for the end of the iterations. When no more iteration can provide an improvement of Xb larger than impmin, the algorithm stops. If not speficief, impmin=1/100.
}
  \item{s}{The threshold used in the iterations. If not specified, the theoretical value of Alquier (2010) is used: s = sqrt(2*v*log(p*K/epsilon)).
}
  \item{v}{The variance of e, if it is known. If not specified, estimated on the data (by a MA(10)-smoothing).
}
}
\value{
\item{Esti.CopyN}{data vector reporting the estimated DNA copy numbers for seleted genes/clones
of all the samples.}
\item{CGH.Array}{a copy of the input data CGH.Array.}
\item{chromosome}{a copy of the input chromosome.}
\item{nucleotide.position}{of copy of the input nucleotide.position.}
\item{FDR}{NULL, FDR is not computed.}
}
\references{
P. Alquier, An Algorithm for Iterative Selection of Blocks of Features, Proceedings of ALT'10, 2010, M. Hutter, F. Stephan, V. Vovk and T. Zeugmann Eds., Lecture Notes in Artificial Intelligence, pp. 35-49, Springer.

R. Tibshirani and P. Wang, Spatial Smoothing and Hot Spot Detection for CGH Data Using the Fused Lasso, Biostatistics, 2008, vol. 9, no. 1, pp 18-29.
}
\author{
Pierre Alquier <alquier@ensae.fr>
}

\examples{
data(CGHDisease1)
cgh = cghISBF(CGHDisease1$data,CGHDisease1$chromosome,CGHDisease1$nucposi,s=1,K=100)
}

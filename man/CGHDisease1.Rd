\name{CGHDisease1}
\alias{CGHDisease1}
\docType{data}
\title{
An example of CGH array.
}
\description{
An example of CGH array for the package ISBF. Taken from the package cghFLasso.
}
\usage{data(CGHDisease1)}
\format{
A list of 3 (data, chr, nucposi).
}
\details{
The package cghFLasso contains several CGH arrays, in two categories (normal,
and disease). This is just a copy of the first array in the disease category.
Quoting the package cghFLasso: the value of each entry is the log fluorescence
ratio resulted from the CGH experiment. The order of the genes/clones in the
rows is the same as the order of the genes/clones on the genome. chr and nucposi
provide chromosome number and nucleotide position for each gene/clone.
}
\source{
http://www-stat.stanford.edu/~tibs/software.html
}
\references{
R. Tibshirani and P. Wang, Spatial Smoothing and Hot Spot Detection for CGH Data Using the Fused Lasso, Biostatistics, 2008, vol. 9, no. 1, pp 18-29.
}
\examples{
data(CGHDisease1)
cgh = cghISBF(CGHDisease1$data,CGHDisease1$chromosome,CGHDisease1$nucposi,s=1,K=100)
}
\keyword{datasets}

\name{isbf}
\alias{isbf}
\title{
Iterative Selection of Blocks of Features - isbf
}
\description{
isbfperforms estimation in the model Y = b + e where the unknown parameter b is sparse, or sparse and constant by blocks. Y is a vector of size p, b a vector of size p and e is the noise.

When b is sparse and constant by blocks, one can use isbf(Y,K=...) where K is the expected maximal size for a block. The method used is Iterative Selection of Blocks of Features procedure of Alquier (2010). When b is only sparse, one can use isbf(Y), as the default value for K is 1. Of course, one can always set K=p, but be careful, the computation time and the memory used is directly proportional to p*K.

NOTE: one can used isbfReg(X,Y) with X the identity matrix instead, but isbf(Y) is really faster.
}
\usage{
isbf(Y, epsilon = 0.05, K = 1, impmin = 1/100,
s = NULL, v = NULL)
}
\arguments{
  \item{Y}{The data. A vector of size p.
}
  \item{epsilon}{The confidence level. The theoretical guarantees in Alquier (2010) is that each iteration of the ISBF procedure gets closer to the real parameter b with probability at least 1-epsilon. When epsilon is very small, the procedure becomes very conservative. When epsilon is too large, there is a risk of overfitting. If not specified, epsilon = 5\%.
}
  \item{K}{The maximal length of blocks. If not specified, K=1, this means we seek for a sparse (not constant by block) parameter b. One should take a larger K is b is really expected to be constant by blocks. If p is quite small (up to 1000), K=p is a reasonnable choice. For larger values of p, please take into account that the computation time and the memory used is directly proportional to p*K.
}
  \item{impmin}{Criterion for the end of the iterations. When no more iteration can provide an improvement of Xb larger than impmin, the algorithm stops. If not speficied, impmin=1/100.
}
  \item{s}{The threshold used in the iterations. If not specified, the theoretical value of Alquier (2010) is used: s = sqrt(2*v*log(p*K/epsilon)).
}
  \item{v}{The variance of e, if it is known. If not specified, estimated on the data (by a MA(10)-smoothing).
}
}
\value{
\item{beta}{The estimated parameter b.}
\item{s}{The value of s.}
\item{impmin}{The value of impmin.}
\item{K}{The value of K.}
}
\references{
P. Alquier, An Algorithm for Iterative Selection of Blocks of Features, Proceedings of ALT'10, 2010, M. Hutter, F. Stephan, V. Vovk and T. Zeugmann Eds., Lecture Notes in Artificial Intelligence, pp. 35-49, Springer.
}
\author{
Pierre Alquier <alquier@ensae.fr>
}

\examples{
# generating data
b = c(rep(0,100),rep(2,40),rep(0,60))
e = rnorm(200,0,0.3)
Y = b + e

# call of isbf
A = isbf(Y,K=200,v=0.3)

# visualization of the results
plot(Y)
lines(A$beta,col="red")
}

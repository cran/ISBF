\name{isbfReg}
\alias{isbfReg}
\title{
Iterative Selection of Blocks of Features in Regression Estimation - isbfReg
}
\description{
isbfReg performs regression estimation in the model Y = Xb + e where the unknown parameter b is sparse, or sparse and constant by blocks. Y is a vector of size n, X a (n,p) matrix, b a vector of size p and e is the noise.

When b is sparse, one can basically use isbfReg(X,Y), the method used is the Iterative Feature Selection procedure of Alquier (2008). When b is sparse and constant by blocks, one can use isbfReg(X,Y,K=...) where K is the expected maximal size for a block. The method used is Iterative Selection of Blocks of Features procedure of Alquier (2010). Of course, one can always set K=p, but be careful, the computation time and the memory used is directly proportional to p*K.
}
\usage{
isbfReg(X, Y, epsilon = 0.05, K = 1, impmin = 1/100, favgroups = 0,
centX = TRUE, centY = TRUE, s = NULL, v = NULL)
}
\arguments{
  \item{X}{The data: the matrix of inputs. Size (n,p).
}
  \item{Y}{The data: the vector of outputs. Size n.
}
  \item{epsilon}{The confidence level. The theoretical guarantees in Alquier (2010) is that each iteration of the ISBF procedure gets closer to the real parameter b with probability at least 1-epsilon. When epsilon is very small, the procedure becomes very conservative. When epsilon is too large, there is a risk of overfitting. If not specified, epsilon = 5\%.
}
  \item{K}{The maximal length of blocks checked in the iterations. If not specified, K=1, this means we seek for a sparse (not constant by block) parameter b, as in Alquier (2008). One should take a larger K is b is really expected to be constant by blocks. If p is quite small (up to 1000), K=p is a reasonnable choice. For larger values of p, please take into account that the computation time and the memory used is directly proportional to p*K.
}
  \item{impmin}{Criterion for the end of the iterations. When no more iteration can provide an improvement of Xb larger than impmin, the algorithm stops. If not speficied, impmin=1/100.
}
  \item{favgroups}{In case of noisy input data, one may want to favor larger groups in order to stabilize estimation. By default, this variable is taken to 0, but take it larger for
noisy input data.
}
  \item{centX}{If TRUE, the function centers the variables in X before processing.
}
  \item{centY}{If TRUE, the function centers the variable Y before processing.
}
  \item{s}{The threshold used in the iterations. If not specified, the theoretical value of Alquier (2010) is used: s = sqrt(2*v*log(p*K/epsilon)).
}
  \item{v}{The variance of e, if it is known. If not specified, this parameter is VERY roughly estimated by var(Y)/2.
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

P. Alquier, Iterative Feature Selection in Least Square Regression Estimation, Annales de l'IHP, B (Proba. Stat.), 2008, vol. 44, no. 1, pp 47-88.
}
\author{
Pierre Alquier <alquier@ensae.fr>
}

\examples{
# generating data
X = matrix(data=rnorm(5000),nr=50,nc=100)
b = c(rep(0,50),rep(-3,30),rep(0,20))
e = rnorm(50,0,0.3)
Y = X\%*\%b + e

# call of isbfReg
A = isbfReg(X,Y,K=100,v=0.3)

# visualization of the results
plot(b)
lines(A$beta,col="red")
}

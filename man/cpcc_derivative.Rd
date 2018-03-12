\name{cpcc_derivative}
\alias{cpcc_derivative}

\title{
Calculates the CPCC derivative for the given input.
}
\description{
Calculates the CPCC derivative for the given input. Should not be used by the user.
}
\usage{
cpcc_derivative(x, var_list)
}
\arguments{
  \item{x}{
    The weight vector of length \code{n-1}, where \code{n} is the number of variables, for the first \code{n-1} . See the reference in \code{\link{find_weights}} for the explanation.
}
  \item{var_list}{
    A list of variables containing the following: data, the relevant data. method, the clustering method. bounds, the bounds on the weights. MMM, Minimal memory mode: see \code{\link{find_weights}} for details. fk, which only exists if MMM is TRUE: see \code{\link{find_weights}} for details.
}
}

\value{
\item{CPCC derivative}{The CPCC derivative of this clustering with the given weights. A vector of length \code{n-1}.}
}
\references{
    Clustering with optimised weights for Gower's metric: Using hierarchical clustering and Quasi-Newton methods to maximise the cophenetic correlation coefficient, Jeroen van den Hoven.
}
\author{
  Jeroen van den Hoven
}

\examples{
data(mtcars)
L = list(data = mtcars, method = "average", bounds = c(0,1), MMM = TRUE)
x = rep(1/11,10)
cpcc_derivative(x,L)
}

\keyword{ cluster }
\keyword{ optimize }

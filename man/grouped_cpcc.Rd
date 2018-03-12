\name{grouped_cpcc}
\alias{grouped_cpcc}

\title{
Calculates the CPCC for the given input with grouped variables.
}
\description{
Calculates the CPCC for the given input with grouped variables. Should not be used by the user.
}
\usage{
grouped_cpcc(x, var_list)
}
\arguments{
  \item{x}{
    The weight vector of length \code{n-1}, where \code{n} is the number of variables, for the first \code{n-1} . See the reference in \code{\link{find_grouped_weights}} for the explanation.
}
  \item{var_list}{
    A list of variables containing the following: data, the relevant data. method, the clustering method. bounds, the bounds on the weights. MMM, Minimal memory mode: see \code{\link{find_grouped_weights}} for details. fk, which only exists if MMM is TRUE: see \code{\link{find_grouped_weights}} for details.
}
}

\value{
\item{CPCC}{The cophenetic correlation coefficient of this clustering with the given weights.}
}
\references{
    Clustering with optimised weights for Gower's metric: Using hierarchical clustering and Quasi-Newton methods to maximise the cophenetic correlation coefficient, Jeroen van den Hoven.
}
\author{
  Jeroen van den Hoven
}

\examples{
data(mtcars)
L = list(data = mtcars, method = "average",
    bounds = c(0,1), MMM = FALSE, combined_indices = as.list(1:ncol(mtcars)))
x = rep(1/11,10)
grouped_cpcc(x,L)
}

\keyword{ cluster }
\keyword{ optimize }

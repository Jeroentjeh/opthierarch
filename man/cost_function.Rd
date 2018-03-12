\name{cost_function}
\alias{cost_function}

\title{
Calculates the cost function of the K-medoids algorithm for the given input.
}
\description{
Calculates the cost function of the K-medoids algorithm for the given input.
 Should not be used by the user.
}
\usage{
cost_function(w, var_list)
}
\arguments{
  \item{w}{
    The weight vector of length \code{n-1}, where \code{n} is the number of variables, for the first \code{n-1} . See the reference in \code{\link{find_kmedoids_weights}} for the explanation.
}
  \item{var_list}{
    A list of variables containing the following: data, the relevant data. method, the number of cluster k. bounds, the bounds on the weights. fk: see \code{\link{find_kmedoids_weights}} for details.
}
}

\value{
\item{Cost function}{The value of the cost function of this clustering with the given weights.}
}
\references{
    Clustering with optimised weights for Gower's metric: Using hierarchical clustering and Quasi-Newton methods to maximise the cophenetic correlation coefficient, Jeroen van den Hoven.
}
\author{
  Jeroen van den Hoven
}

\examples{
data(faithful)
L = list(data = faithful, k = 2, bounds = c(0,1))
w = 0.4
cost_function(w,L)
}

\keyword{ cluster }
\keyword{ optimize }

\name{cost_function_derivative}
\alias{cost_function_derivative}

\title{
Calculates the derivative for the cost function of the K-medoids algorithm for the given input.
}
\description{
Calculates the derivative for the cost function of the K-medoids algorithm for the given input. Should not be used by the user.
}
\usage{
cost_function_derivative(w, var_list)
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
\item{Cost function derivative}{The cost function derivative of this clustering with the given weights. A vector of length \code{n-1}.}
}
\references{
    Clustering with optimised weights for Gower's metric: Using hierarchical clustering and Quasi-Newton methods to maximise the cophenetic correlation coefficient, Jeroen van den Hoven.
}
\author{
  Jeroen van den Hoven
}

\keyword{ cluster }
\keyword{ optimize }

\name{opt_kmedoids}
\alias{opt_kmedoids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find optimal weights for Gower's metric and return the kmedoids clustering
}
\description{
This function is a wrapper function for find_kmedoids_weights that returns both the kmedoids clustering and the result from find_kmedoids_weights.
}
\usage{
opt_kmedoids(data, k, start_values = rep(1 / ncol(data), ncol(data) - 1),
                        n_iterate = 10,
                        bounds=c(1 / (3 * ncol(data)),
                        1 - (ncol(data) - 1) / (3 * ncol(data))))
                        }
\arguments{
  \item{data}{the data that needs to be clustered, provided as a dataframe or (numeric) matrix. It is assumed that rows correspond to instances and columns correspond to features.}
  \item{start_values}{a vector containing the initial values of the weights. Defaults to \code{1/ncol(data)} for all variables. They must not be negative. Furthermore, \code{1 - sum(start_values) >= bounds[1]} must hold. Due to the way the algorithm is programmed, you only have to supply values for the first ncol(data) - 1 variables!}
  \item{k}{the number of clusters, the K in K-medoids.}
  \item{n_iterate}{the maximum number of iterations used by the quasi-newton method \code{optim}. Defaults to 10.}
  \item{bounds}{a vector of size 2 containing the lower and upper bound in position 1 and 2 respectively. The lower bound must not be lower than 0 and not higher than \code{1/ncol(data)}.
The upper bound will be set to the minimum of its current value and \code{1-(ncol(data)-1)*bounds[1]}. For more information, see details.}
}

\details{
See \code{\link{find_kmedoids_weights}} for more details.
}
\value{
This function returns a list of two named lists:
\item{clustering}{the kmedoids clustering. See \code{\link{pam}} for more information.}
\item{opt_result}{the results of the optimisation algorithm. See \code{\link{find_kmedoids_weights}} for more information.}

}
\references{
    Clustering with optimised weights for Gower's metric: Using hierarchical clustering and Quasi-Newton methods to maximise the cophenetic correlation coefficient, Jeroen van den Hoven.
}
\author{
  Jeroen van den Hoven
}
\note{
This package requires the cluster package.
}

\seealso{
\code{\link{find_kmedoids_weights}}
}
\examples{
## Basic example
data(mtcars)
find_kmedoids_weights(mtcars, k = 3)

## Using custom bounds
find_kmedoids_weights(mtcars, k = 3, bounds = c(0.05, 0.2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }

\name{find_grouped_weights}
\alias{find_grouped_weights}
\title{Optimise Gower's weights for hierarchical clustering}

\description{
  Uses the limited memory BFGS algorithm with bounds from \code{optim} to optimise the weights of a hierarchical clustering with the goal of maximising the cophenetic correlation coefficient. This function gives you the chance
  to set groups of variables that should retain the same weight.
}
\usage{
find_grouped_weights(data, combined_indices = as.list(1:ncol(data)),
                                start_values = rep(1 / ncol(data),
                                length(combined_indices) - 1),
                                n_iterate = 10, clust_method = "average",
                                bounds = c(1 / (3 * ncol(data)),
                                1 - (ncol(data) - 1) / (3 * ncol(data))),
                                minimal_memory_mode = F)}

\arguments{
  \item{data}{the data that needs to be clustered, provided as a dataframe or (numeric) matrix. It is assumed that rows correspond to instances and columns correspond to features.}
  \item{combined_indices}{the indices that are combined into groups. Present this as a list. Each list entry should contain, in vector format, the indices of the variables that are contained in that group}
  \item{start_values}{a vector containing the initial values of the weights. Defaults to \code{1/ncol(data)} for all variables. They must not be negative. Furthermore, \code{1 - sum(start_values) >= bounds[1]} must hold. Due to the way the algorithm is programmed, you only have to supply values for the first ncol(data) - 1 variables! }
  \item{n_iterate}{the maximum number of iterations used by the quasi-newton method \code{optim}. Defaults to 10.}
  \item{clust_method}{a string containing the type of linkage function used by \code{hclust}. Defaults to average linkage}
  \item{bounds}{a vector of size 2 containing the lower and upper bound in position 1 and 2 respectively. The lower bound must not be lower than 0 and not higher than \code{1/ncol(data)}.
The upper bound will be set to the minimum of its current value and \code{1-(ncol(data)-1)*bounds[1]}. For more information, see details.}
  \item{minimal_memory_mode}{logical that determines whether the algorithm calculates the differences for each instance and each variables beforehand or calculates them live each time. The first will be chosen when this variable is FALSE,
the second one will be chosen when this variable is TRUE. Note that this requires k vectors of size \code{n(n-1)/2} to be stored, where k is the number of columns and n the number of rows. If you have the required memory, setting this
to FALSE is definitely worth the speed increase, which seems to be a factor of something between 2 and 3.}
}

\details{
Contrary to intuition, the start_values vector should not be equal to the number of columns in data. It should have one weights less than the number of groups. The reason for this is that we set the sum of all weights to equal a constant (1 in this case), allowing us to not have to set a weight for one group. This allows us to skip the calculations for that variable, saving some time. The weight for this variable should of course still abide the given bounds!
}
\value{
The result is the output of \code{\link{optim}}.

\item{par}{The best set of parameters found.}
\item{value}{The value of fn corresponding to par.}
\item{counts}{A two-element integer vector giving the number of calls to fn and gr respectively. This excludes those calls needed to compute the Hessian, if requested, and any calls to fn to compute a finite-difference approximation to the gradient.}
\item{convergence}{An integer code. 0 indicates successful completion (which is always the case for "SANN" and "Brent"). Possible  error codes are noted in the \code{\link{optim}} documentation}
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


\examples{
## Basic example
data(esoph)
find_grouped_weights(esoph)

## Using custom groups of weights
find_grouped_weights(esoph, list(c(1,5),2,3,4))
}

\keyword{ cluster }
\keyword{ optimize }

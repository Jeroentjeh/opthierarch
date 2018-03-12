\name{opt_grouped_hierarchical}
\alias{opt_grouped_hierarchical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find optimal weights for Gower's metric and return the corresponding clustering
}
\description{
This function is a wrapper function for find_grouped_weights that returns both the hierarchical clustering and the result from find_grouped_weights.
}
\usage{
opt_grouped_hierarchical(data, combined_indices = as.list(1:ncol(data)),
            start_values = rep(1 / ncol(data), length(combined_indices) - 1),
            n_iterate = 10, clust_method = "average",
            bounds = c(1 / (3 * ncol(data)),1 - (ncol(data) - 1) / (3 * ncol(data))),
            minimal_memory_mode = F)}
\arguments{
The arguments are the same as those of \code{\link{find_grouped_weights}}: see that for details.

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
See \code{\link{find_grouped_weights}} for more details.
}
\value{
This function returns a list of two named lists:
\item{clustering}{the hierarchical clustering. See \code{\link{hclust}} for more information.}
\item{opt_result}{the results of the optimisation algorithm. See \code{\link{find_grouped_weights}} for more information.}

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
\code{\link{find_grouped_weights}}
}
\examples{
## Basic example
data(mtcars)
opt_grouped_hierarchical(mtcars)

## Using custom bounds
opt_grouped_hierarchical(mtcars, combined_indices = (c(list(c(1,2)),as.list(3:11))),
        bounds = c(0, 1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_homo.R
\name{is_homo}
\alias{is_homo}
\title{U-statistic based homogeneity test}
\usage{
is_homo(md = NULL, data = NULL, rep = 10)
}
\arguments{
\item{md}{Matrix of distances between all data points.}

\item{data}{Data matrix. Each row represents an observation.}

\item{rep}{Number of times to repeat optimization procedure. Important for problems with
multiple optima.}
}
\value{
Returns a list with the following elements:\describe{
  \item{minFobj}{Test statistic. Minimum of the objective function for optimization (-stdBn).}
  \item{group1}{Elements in group 1 in the maximal partition. (obs: this is not the best
  partition for the data, see \code{uclust})}
  \item{group2}{Elements in group 2 in the maximal partition.}
  \item{p.MaxTest}{P-value for the homogeneity test.}
  \item{Rep.Fobj}{Values for the minimum objective function on all \code{rep} optimization runs.}
  \item{bootB}{Resampling variance estimate for partitions with groups of size n/2 (or (n-1)/2 and (n+1)/2 if n is odd).}
  \item{bootB1}{Resampling variance estimate for partitions with one group of size 1.}

}
}
\description{
Homogeneity test based on the statistic \code{bn}. The test assesses whether there exists a data partition
for which group separation is statistically significant according to the U-test. The null hypothesis
is overall sample homogeneity, and a sample is considered homogeneous if it cannot be divided into
two statistically significant subgroups.
}
\details{
This is the homogeneity test of Cybis et al. (2017) extended to account for groups of size 1.
The test is performed through two steps: an optimization procedure that finds the data partition that
maximizes the standardized Bn and a test for the resulting maximal partition. Should be used in high dimension small sample size settings.



Either \code{data} or \code{md} should be provided.
If data are entered directly, Bn will be computed considering the squared Euclidean distance.

Variance of \code{bn} is estimated through resampling, and thus, p-values may vary a bit in different runs.

For more detail see Cybis, Gabriela B., Marcio Valk, and Sílvia RC Lopes. "Clustering and classification problems in genetics through U-statistics."
Journal of Statistical Computation and Simulation 88.10 (2018)
and Valk, Marcio, and Gabriela Bettella Cybis. "U-statistical inference for hierarchical clustering." Journal of Computational and Graphical Statistics 30(1) (2021).
}
\examples{
x = matrix(rnorm(500000),nrow=50)  #creating homogeneous Gaussian dataset
res = is_homo(data=x)

x[1:30,] = x[1:30,]+0.15   #Heterogeneous dataset (first 30 samples have different mean)
res = is_homo(data=x)

md = as.matrix(dist(x)^2)   #squared Euclidean distances for the same data
res = is_homo(md)

# Multidimensional sacling plot of distance matrix
fit <- cmdscale(md, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
plot(x,y, main=paste("Homogeneity test: p-value =",res$p.MaxTest))

}

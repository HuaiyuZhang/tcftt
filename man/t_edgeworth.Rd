% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{t_edgeworth}
\alias{t_edgeworth}
\title{Edgeworth expansion for Welch's t-statistic}
\usage{
t_edgeworth(
  x,
  order = 2,
  n1,
  n2,
  mu1,
  mu2,
  sigma1,
  sigma2,
  gamma1,
  gamma2,
  tau1,
  tau2
)
}
\arguments{
\item{x}{a real number.}

\item{order}{the order of edgeworth expansion. Valid options are 0, 1, and 2. If set to 0,
it reduces to approximation based on the central limit theorem and returns the CDF of standard normal distribution
evaluated at x.}

\item{n1}{sample size for the sample from the first population.}

\item{n2}{sample size for the sample from the second population.}

\item{mu1}{mean of the first population.}

\item{mu2}{mean of the second population.}

\item{sigma1}{standard deviation of the first population.}

\item{sigma2}{standard deviation of the second population.}

\item{gamma1}{skewness of the first population.}

\item{gamma2}{skewness of the second population.}

\item{tau1}{kurtosis of the first population.}

\item{tau2}{kurtosis of the second population.}
}
\value{
Edgeworth expansion evaluated at x.
}
\description{
This function provides approximation for the cumulative distribution function
of the sampling distribution of the Welch's t-statistic using Normal distribution, first order or second order Edgeworth expansion.
}
\examples{
t_edgeworth(1.96, order=2,
n1=20, n2=30,
mu1=0, mu2=0,
sigma1=1, sigma2=0.5,
gamma1=1, gamma2=0,
tau1=6, tau2=0)
}

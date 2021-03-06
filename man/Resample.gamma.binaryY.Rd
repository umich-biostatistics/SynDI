% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ResampleGamma.R
\name{Resample.gamma.binaryY}
\alias{Resample.gamma.binaryY}
\title{Resample for bootstrap variance for binary Y}
\usage{
Resample.gamma.binaryY(data, indices)
}
\arguments{
\item{data}{synthetic data}

\item{indices}{row indices to replicate}
}
\value{
numeric vector of regression coefficients
}
\description{
Resampling function to get bootstrap variance for binary Y. Note that
readers need to modify the existing function Resample.gamma.binaryY() to
match their own Steps 1-5. It was only included in the package for the
purpose of providing an example.
}
\references{
Reference: Gu, T., Taylor, J.M.G. and Mukherjee, B. (2021) Regression
inference for multiple populations by integrating summary-level data using stacked
imputations \url{https://arxiv.org/abs/2106.06835}.
}

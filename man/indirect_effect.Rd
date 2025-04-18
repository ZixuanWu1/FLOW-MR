% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other_effects.R
\name{indirect_effect}
\alias{indirect_effect}
\title{Indirrect effedt}
\usage{
indirect_effect(
  results,
  K,
  warmup = 3000,
  path = NULL,
  quantiles = c(0.025, 0.5, 0.975)
)
}
\arguments{
\item{results}{results from gibbs_wrapper}

\item{K}{number of phenotypes}

\item{quantiles}{A vector of quantiles of interests}

\item{warump}{The length of warmup period. Default is 3000}
}
\value{
A K by K by 4 array of summary of direct effect.
}
\description{
Estimate the indirrect effect of an exposure on the outcome
}

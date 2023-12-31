% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VisualFuns.R
\name{clusGapDiscrMDS}
\alias{clusGapDiscrMDS}
\title{Multidimensional scaling plot for categorical data.
Scatter dot dimension reduction plot obtained from user input number of clusters `nCl`, specified categorical distance `distName` and clustering method function `FUNcluster.
First two dimensions are used including the estimated percentage explained variability. Dots are color-coded according to the assigned cluster.
The function outputs a ggplot object where additional layers can be added for customization purposes.}
\usage{
clusGapDiscrMDS(x, distName, nCl, FUNcluster, title = NULL)
}
\arguments{
\item{x}{matrix object or data.frame}

\item{distName}{Name of categorical distance to apply.
Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'.}

\item{nCl}{number of clusters to plot}

\item{FUNcluster}{a function that accepts as first argument a matrix like `x`; second argument specifies number of `k` (k=>2) clusters
This function returns a list with a component named `cluster`, a vector of length `n=nrow(x)` of integers from `1:k` indicating observation cluster assignment.}

\item{title}{Plot title}
}
\value{
ggplot object.
}
\description{
Multidimensional scaling plot for categorical data.
Scatter dot dimension reduction plot obtained from user input number of clusters `nCl`, specified categorical distance `distName` and clustering method function `FUNcluster.
First two dimensions are used including the estimated percentage explained variability. Dots are color-coded according to the assigned cluster.
The function outputs a ggplot object where additional layers can be added for customization purposes.
}

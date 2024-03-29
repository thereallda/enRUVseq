% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{ggDotPlot}
\alias{ggDotPlot}
\title{Dot-plot with mean_sd bar}
\usage{
ggDotPlot(data, x, y, fill, palette = NULL)
}
\arguments{
\item{data}{A data.frame (or a tibble).}

\item{x}{The grouping variable from the \code{data}.}

\item{y}{The value variable from the \code{data}.}

\item{fill}{The fill variable from the \code{data}.}

\item{palette}{The fill palette for different groups.}
}
\value{
ggplot2 object
}
\description{
Dot-plot with mean_sd bar
}

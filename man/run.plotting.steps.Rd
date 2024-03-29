% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pause_workflow_funcs_v2.R
\name{run.plotting.steps}
\alias{run.plotting.steps}
\title{run the plotting steps given a pause.body object (with desired genes) and pair of bigWigs from two conditions being compared}
\usage{
run.plotting.steps(
  df.pausebody,
  cond1.name,
  bw.cond1.plus,
  bw.cond1.minus,
  cond2.name,
  bw.cond2.plus,
  bw.cond2.minus,
  factor.name = "efferocytosis",
  color.names = c(rgb(0, 0, 1, 1/2), rgb(1, 0, 0, 1/2))
)
}
\arguments{
\item{df.pausebody}{an object containing start, end coordinates of pause window and gene body.}

\item{cond1.name}{name of baseline condition}

\item{bw.cond1.plus}{bigWig file containing plus strand data from baseline condition}

\item{cond2.name}{name of the condition which is being compared against the baseline condition}

\item{bw.cond2.plus}{bigWig file containing plus strand data from the condition which is being compared against baseline condition}

\item{bw.cond2.minus}{bigWig file containing minus strand data from the condition which is being compared against baseline condition'}

\item{factor.name}{descriptive name of the factor involved in the comparison of condition to baseline}

\item{color.names}{to be used as colors for plotting. The order of colors chosen will be same as the alphabetical order of names of two conditions}

\item{bw, cond1.minus}{bigWig file containing minus strand data from baseline condition}
}
\description{
run the plotting steps given a pause.body object (with desired genes) and pair of bigWigs from two conditions being compared
}

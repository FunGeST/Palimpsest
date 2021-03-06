% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/palimpsest_input.R
\name{makeMutypeMatFromVcf}
\alias{makeMutypeMatFromVcf}
\title{makeMutypeMatFromVcf}
\usage{
makeMutypeMatFromVcf(
  vcf,
  sample.col = "Sample",
  mutcat.col = "SBS_cat3",
  mutypes = c("CA", "CG", "CT", "TA", "TC", "TG"),
  proportion = TRUE
)
}
\arguments{
\item{vcf}{vcf data frame containing the mutations/SVs}

\item{sample.col}{sample.col}

\item{mutcat.col}{mutcat.col}

\item{mutypes}{mutypes}

\item{proportion}{proportion}
}
\description{
Function to create a matrix in mutation type x sample format with either counts or proportions
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FREEC.R
\name{parseCfreec_file}
\alias{parseCfreec_file}
\title{Load a ControlFreec _CNV output from a tab delimited file
this function can also use a _CNV file with pvalues added by the accompanying script in the controlFreec package.
it expects to find a tab delimited file
The columns are:
\itemize{
\item chromosome
\item start position
\item end position
\item copy.number
\item status
\item genotype
\item uncertainty of the assignment
\item WilcoxonRankSumTestPvalue
\item KolmogorovSmirnovPvalue
Note: at least for v 11.5 it is not possible to understand how subclones are coded in the output files.
As a result we only parse the main clone information (i.e no subclones). A related ticket has been
opened on GitHub to help understand the C-FREEC outputs (https://github.com/BoevaLab/FREEC/issues/84)
}}
\usage{
parseCfreec_file(cfreec.fn)
}
\arguments{
\item{cfreec.fn}{input filename with CNV segments,}

\item{filterLowUncerntainty}{if set to TRUE segments with uncertainty equal to -1 will be discarded}

\item{pvalThreshold}{threshold of pvalues to filter if the pvalues are included in the file}
}
\value{
A GRanges object with the segments predicted by ControlFreec
}
\description{
Load a ControlFreec _CNV output from a tab delimited file
this function can also use a _CNV file with pvalues added by the accompanying script in the controlFreec package.
it expects to find a tab delimited file
The columns are:
\itemize{
\item chromosome
\item start position
\item end position
\item copy.number
\item status
\item genotype
\item uncertainty of the assignment
\item WilcoxonRankSumTestPvalue
\item KolmogorovSmirnovPvalue
Note: at least for v 11.5 it is not possible to understand how subclones are coded in the output files.
As a result we only parse the main clone information (i.e no subclones). A related ticket has been
opened on GitHub to help understand the C-FREEC outputs (https://github.com/BoevaLab/FREEC/issues/84)
}
}

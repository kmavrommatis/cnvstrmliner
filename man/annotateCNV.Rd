% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{annotateCNV}
\alias{annotateCNV}
\title{Annotate a CNV}
\usage{
annotateCNV(major.cn = NA, minor.cn = 0, total.cn = NA, na.minor = 0)
}
\arguments{
\item{major.cn}{Major Allele copy number (integer)}

\item{minor.cn}{Minor Allele copy number (integer)}

\item{total.cn}{Total copy number}

\item{na.minor}{The default value that the minor allele will take if it is NA. Defaults to 0}
}
\value{
Returns a string indicating the type of the copy number aberration.
}
\description{
annotate a CNV segment with a flag describing a Copy number aberration according to the table
}
\details{
\tabular{llll}{
   total cn \tab minor cn \tab type \tab meaning of type \cr
   0 \tab 0 \tab HOMDEL \tab Homozygous deletion \cr
   0 \tab NA \tab HOMDEL \tab Homozygous deletion \cr
   1 \tab 0 \tab HETDEL \tab Heterozygous deletion (1 out of 2 alleles) \cr
   1 \tab NA \tab HETDEL \tab Heterozygous deletion (this is an edge case that should not occur). \cr
   2 \tab 0 \tab NEUT-LOH \tab Loss of heterozygosity without change in the number of alleles \cr
   2 \tab 1 \tab NEUT \tab Normal \cr
   2 \tab NA \tab NEUT/Unknown \tab Normal or undetermined by the algorithm \cr
   2+ \tab 0 \tab AMP-LOH \tab Amplification with loss of heterozygosity \cr
   2+ \tab 1+ \tab AMP \tab Amplification \cr
   2+ \tab NA \tab AMP/LOH. \tab Amplification with (possible) loss of heterozygosity. Cannot determine. \cr
}
}
\section{Notes}{

The major.cn and minor.cn should add up to total.cn.
Either total.cn or major.cn are necessary.
If one of these is not provided the function
will calculate it assuming total.cn=major.cn + minor.cn
Minor.cn is not necessary, it defaults to 0 meaning that there is no minor allele
}

\section{Ploidy consideration}{

This function assumes that normal ploidy is 2.
}


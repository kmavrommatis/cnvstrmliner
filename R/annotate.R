#' Annotate a CNV
#'
#' annotate a CNV segment with a flag describing a Copy number aberration according to the table
#'
#' |total cn | minor cn |	type         | meaning of type                                                        |
#' |---------|----------|--------------|------------------------------------------------------------------------|
#' |0        |	0 	    | HOMDEL       | Homozygous deletion                                                    |
#' |0        |	NA 	    | HOMDEL       | Homozygous deletion                                                    |
#' |1        |	0 	    | HETDEL       | Heterozygous deletion (1 out of 2 alleles)                             |
#' |1        |	NA 	    | HETDEL       | Heterozygous deletion (this is an edge case that should not occur).    |
#' |2        |	0 	    | NEUT-LOH     | Loss of heterozygosity without change in the number of alleles         |
#' |2        |	1 	    | NEUT         | Normal                                                                 |
#' |2        |	NA 	    | NEUT/Unknown | Normal or undetermined by the algorithm                                |
#' |2+       |	0 	    | AMP-LOH      | Amplification with loss of heterozygosity                              |
#' |2+       |	1+ 	    | AMP          | Amplification                                                          |
#' |2+       |	NA 	    | AMP/LOH.     | Amplification with (possible) loss of heterozygosity. Cannot determine.|
#'
#'
#' @param major.cn Major Allele copy number (integer)
#' @param minor.cn Minor Allele copy number (integer)
#' @param total.cn Total copy number
#' @param na.minor  The default value that the minor allele will take if it is NA. Defaults to 0
#' @section Notes:
#' The major.cn and minor.cn should add up to total.cn.
#' Either total.cn or major.cn are necessary.
#' If one of these is not provided the function
#' will calculate it assuming total.cn=major.cn + minor.cn
#' Minor.cn is not necessary, it defaults to 0 meaning that there is no minor allele
#' @section Ploidy consideration:
#' This function assumes that normal ploidy is 2.
#' @return Returns a string indicating the type of the copy number aberration.
#'
#' @export
annotateCNV=function(major.cn=NA, minor.cn=0, total.cn=NA, na.minor=0){
  #message(major.cn, " ", class(major.cn)," " ,minor.cn," ",class(minor.cn), " " , total.cn)
  if(is.na( total.cn ) & is.na(major.cn)){
    stop("Please provide either the total number of alleles or the major allele copies")
  }
  major.cn=as.numeric(major.cn)
  minor.cn=as.numeric(minor.cn)
  total.cn=as.numeric(total.cn)
  if(is.na(total.cn)){
    total.cn=ifelse( is.na(minor.cn), major.cn, major.cn + minor.cn)
  }
  if(total.cn <0){
    stop("total.cn (major.cn + minor.cn) cannot be a negative number")
  }
  if(is.na(minor.cn)){
    #message("The minor.cn is not defined. This is a tricky situation in most cases, we default to 0 minor allele copies, but can be changed with the na.minor argument")
    minor.cn=na.minor
  }
  if(minor.cn <0){
    stop("minor.cn (major.cn + minor.cn) cannot be a negative number")
  }

  flag="UNKNOWN"
  # set the flag
  if( total.cn == 0 & (minor.cn ==0 | is.na(minor.cn))){ flag="HOMDEL"}
  if( total.cn == 1 & (minor.cn ==0 | is.na(minor.cn))){ flag="HETDEL"}
  if( total.cn == 2 & (!is.na(minor.cn) & minor.cn ==0 )){ flag="NEUT-LOH"}
  if( total.cn == 2 & (!is.na(minor.cn) & minor.cn ==1 )){ flag="NEUT"}
  if( total.cn == 2 & is.na(minor.cn) ){ flag="NEUT/UNKNOWN"}
  if( total.cn >  2 & (!is.na(minor.cn) & minor.cn ==0 )){ flag="AMP-LOH"}
  if( total.cn >  2 & (!is.na(minor.cn) & minor.cn >= 1 )){ flag="AMP"}
  if( total.cn >  2 & is.na(minor.cn)  ){ flag="AMP/LOH"}

  return( flag )

}

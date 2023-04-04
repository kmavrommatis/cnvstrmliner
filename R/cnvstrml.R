
#' check the validity of the inputs for the cnvstrml object
#'
#' @param object the cnvstrml object ot check
#' @return TRUE or errors
check_cnvstrml = function(object) {
  errors <- character()
  if(object@purity <= 0 ){ msg="Purity cannot be a negative number. It has to be between 0 and 1"; errors=c(errors,msg)}
  if(object@purity >1 ){ msg="Purity cannot be a greater than 1. It has to be between 0 and 1"; errors=c(errors,msg)}
  if(object@ploidy <=0 ){ msg="Ploidy cannot be negative or 0. It has to be a positive number"; errors=c(errors,msg)}
  if(! object@genome %in% c("hg38","hg19") ){msg="The genome should be either hg19 or hg38"; errors=c(errors,msg)}
  if( !is.na(object@sex)){
    object@sex=toupper(object@sex)
    if( ! object@sex %in% c( 'XY','XX','M','F') ){
      msg="Not a valid sex. Sex can be either XY,M for male, or XX, F for female"; errors=c(errors,msg)
    }
  }
  if (length(errors) == 0) TRUE else errors
}



#setOldClass("GenomicRanges::GRanges")
#' create the a CNVstrml object
#'
#' cnvstrml is a structure that holds the information for the CNV predictions
#' and some additional information provided by the methods
#'
#' @slot segments set of CNV segments (GRanges) as predicted by a CNV calling method
#' @slot purity (real) the purity of the sample (default 1)
#' @slot ploidy (real) the ploidy of the sample (default 2)
#' @slot genome (character) the name of the genome to use
#' @slot sample (character) name of sample
#' @slot method (character) the name of the method used to produce the data
#' @slot sex (character) the sex of the sample
#' @slot segtype (character) the type of data stored in this object
#' @docType class
#' @export
#' @section Notes:
#' This is the main object provided by this package.
#'
#' Contains information about the CNV segments
#' the purity and ploidy of the samples,
#' and metadata about the method name, sample name and the
#' reference genome used
#'
cnvstrml=setClass("cnvstrml",
                  #contains=c("GenomicRanges::GRanges"),
                  slots=c(
                          segments="GRanges",
                          purity="numeric",
                          ploidy="numeric",
                          method="character",
                          sample="character",
                          genome="character",
                          sex="character",
                          segtype="character"),
                  prototype=list(
                          purity=1,
                          ploidy=2,
                          method=NA_character_,
                          sample=NA_character_,
                          genome='hg38',
                          sex=NA_character_,
                          segtype='single sample'
                  ),
                  validity = check_cnvstrml
                  )


#' initialize the object
#'
#' check that purity is within limits
#' updates the sex based on the presence of chrY events if the user has not prodived that
#'
#' @param .Object an .Objectect of class cnvstrml
#' @param ... additional arguments
#' @return an initialized .Objectect
#'
setMethod("initialize", "cnvstrml",
          function(.Object, ...) {
            .Object <- callNextMethod()
            #callNextMethod()

            # check and infer the sex
            if(is.na(.Object@sex)){
                sex=inferSex( .Object)
                .Object@sex=toupper(sex)
            }
            .Object
          }
)


#' show the contents of the object
#'
#' @param object Object of type cnvstrml
#'
#'
setMethod("show", "cnvstrml",
          function(object) {
            cat("Contents of sample ", object@sample,"\n")
            cat(paste("Method used to produce data ", object@method ,"\n"))
            cat(paste("Purity is set to ", object@purity, "\n"))
            cat(paste("Ploidy is set to ", object@ploidy, " \n"))
            cat(paste("Sex is set to ", object@sex, "\n"))
            cat(paste(length(object@segments)," segments\n"))
          }
)

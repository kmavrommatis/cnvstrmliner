#' Get the  ploidy for the sample
#'
#' Two parameters can be retrieved:
#'
#' **rounded ploidy** is the rounded ploidy (e.g. 2.1 will become 2). This is the default.
#'
#' **ploidy** is the actual value provided by the user
#'
#'
#' @param .Object object of type cnvstrml,
#' @param rounded get the rounded ploidy (default TRUE)
#' @docType methods
#' @return the ploidy
#'
setGeneric("getPloidy", function(.Object, rounded){
  standardGeneric("getPloidy")
})

#' Get the  ploidy for the sample
#'
#' @rdname getPloidy
#' @name getPloidy
#' @export getPloidy
#'
setMethod("getPloidy",
          signature(.Object="cnvstrml", rounded="logical"),
          function(.Object, rounded){
            if(rounded == FALSE){
              return( .Object@ploidy)
            }
            return(round(.Object@ploidy) )
          }
)

#' Get the  ploidy for the sample
#'
#' @rdname getPloidy
#' @name getPloidy
#' @export getPloidy
#'
setMethod("getPloidy",
          signature(.Object="cnvstrml"),
          function(.Object, rounded){
            return( getPloidy( .Object, TRUE ))
          }
)

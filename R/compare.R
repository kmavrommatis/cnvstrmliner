#' compare the segments between two cnvstrml objects
#'
#' Compares the query to the subject segments
#' Finds teh common ranges, i.e. overlapping using findOverlaps
#' it performs the comparisons at three levels
#' * cnv status (e.g. AMP, NEUT, DEL etc)
#' * cnv copies (integer number of copies). This is applied only when status TP
#' * cnv copies and ccf. CCF is adjusted by wiggle. This is applied only when status is TP
#' TP
#' @param query (cnvstrml)
#' @param subject (cnvstrml)
#' @param wiggle (float) the difference between ccfs that are allowed to be considered the same (e.g. with wiggle =0.1 ccf 0.3 and 0.39 are considered the same)
#' @return a table with numbers of TP, FP, FN
#' @export

compareSegmentStatus=function( query,subject, wiggle=0.1){
  query.gr=query@segments
  subject.gr=subject@segments

  queryOv=pintersect( findOverlapPairs( query.gr, subject.gr ))
  subjectOv=pintersect( findOverlapPairs( subject.gr, query.gr ))

  keep=intersect( which(subjectOv$hit==TRUE),which(queryOv$hit==TRUE) )

  queryOv=queryOv[keep]

  subjectOv=subjectOv[keep]
  # all the ranges
  all=length( subjectOv )
  all_bp=sum(width(subjectOv))

  # find the common ranges (same cnv.flag)
  tt=cbind(
    chrom=seqnames(queryOv)%>% as.character(),
    start=start(queryOv),
    end=end(queryOv),
    query=queryOv$cnv.annotation,
    subject=subjectOv$cnv.annotation,
    query_CNV=queryOv$total.allele.copies,
    subject_CNV=subjectOv$total.allele.copies,
    query_CCF=queryOv$cnv.ccf,
    subject_CCF=subjectOv$cnv.ccf
  ) %>% as.data.frame()
  tt=tt %>%
    dplyr::mutate(across( c(query_CNV,subject_CNV,query_CCF, subject_CCF) , as.numeric ) )  %>%
    dplyr::mutate( status= ifelse( query == subject, "TP", ifelse( query == "NEUT", "FP" ,"FN"))) %>%
    dplyr::mutate( cnv =   ifelse( status == "TP", ifelse( query_CNV == subject_CNV, 'TP', ifelse( query_CNV < subject_CNV, 'FP', 'FN'))  , NA) ) %>%
    dplyr::mutate( ccf =   ifelse( status == "TP", ifelse( abs(query_CCF - subject_CCF)<wiggle, 'TP',ifelse( query_CCF< subject_CCF, 'FP', 'FN'  )), NA)  )

  return(tt)

}

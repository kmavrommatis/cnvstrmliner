
#' get the consensus sequence from a list of cnvstrml objects that represent the CNV segments
#'
#' The consensus is a new Granges that contains only the shared segments
#' The metadata of the consensus have the cnv.annotations for each of the elements of the segments
#' @section Note
#' All the elements of the segmentList should have the same genome and sample name.
#'
#' @param segmentsList (list) a list of cnvstrml objects
#' @param ignore.samples (logical) if set to FALSE then the sample name of each element in segmentsList has to be the same (default FALSE)
#' @param method (character)
#' @returns a cnvstrml with all the segments that are shared across all samples.
#' @export

findConsensus=function(segmentsList,  ignore.samples=FALSE, method=c("qualitative","quantitative")){
  method=match.arg(method)
  samples=sapply( segmentsList, function(x){x@sample}) %>% unique()
  genomes=sapply( segmentsList, function(x){x@genome}) %>% unique()
  sex=sapply( segmentsList, function(x){ x@sex}) %>% unique()
  if('F' %in% sex){ sex = 'F'}else{sex='M'}
  errormesg=c()
  error=FALSE
  if(ignore.samples == FALSE & length(samples)>1){ msg="The samples in this list have different names. All the elements of the list should have the same sample name."; errormesg=c(errormesg, msg); error=TRUE}
  if(length(genomes)>1){ msg="The samples in this list have different genomes All the elements of the list should have the same genome";  errormesg=c(errormesg, msg); error=TRUE}
  if( error==TRUE){
    stop(paste(errormesg, "\n"))
  }

  segsList=lapply( segmentsList , function(x){x@segments})
  consensus=findConsensusGR( segsList , method)

  if(ignore.samples==FALSE & length(names(segmentsList) == length(segmentsList))){
    colnames=names(segmentsList) %>% make.unique()
    colnames(S4Vectors::mcols( consensus ))= c( colnames, 'majority')
  }

  purities=sapply( segmentsList, function(x){x@purity}) %>% mean(. , na.rm=TRUE)
  ploidies=sapply( segmentsList, function(x){x@ploidy}) %>% mean(. , na.rm=TRUE)
  methods=sapply( segmentsList, function(x){x@method}) %>% unique() %>% paste(.,collapse=",")



  obj=cnvstrml(
    segments=consensus,
    purity = purities,
    ploidy = ploidies,
    method= methods,
    sample=samples,
    sex=sex,
    segtype=paste('consensus',method),
    genome=genomes
  )
  obj

}

#' fill the gaps between the segments with segments of neutral copy number
#'
#' This function finds the gaps between the segments that have been called
#' and adds them as new segments with neutral copy number.
#' It is suggested to run this function _after_ applying the consensus with the function applyConsensus()
#'
#' @param segments (cnvstrml) the object which will be modified
#' @param ploidy (integer) expected copy number of the neutral region (default=2)
#' @param major.allele.copies (integer) expected number of the major allele in the neutral regions (default =1)
#' @return an updated copy of the input object
#' @export
fillGaps=function( segments, ploidy=2, major.allele.copies=1){

  segs=segments@segments
  gaps=GenomicRanges::gaps( segs)
  gaps=gaps[ strand(gaps)=="*"]
  if(length(gaps)==0){return( cnv )}

  S4Vectors::mcols( gaps)=S4Vectors::mcols( gaps) %>%data.frame( ) %>%
    dplyr::mutate(total.allele.copies=ploidy,
                  major.allele.copies=major.allele.copies,
                  minor.allele.copies=total.allele.copies-major.allele.copies,
                  cnv.ccf=1,
                  cnv.annotation='NEUT',
                  cnv.loh=FALSE,
                  cnv.lrr=0)

  segments@segments=c(segs, gaps) %>% GenomeInfoDb::sortSeqlevels() %>% sort
  return(segments)
}



#' apply the consensus to the segments
#'
#' this function splits each GRanges object in the segments
#' to the consensus segment and the unique part
#' The 'hit' column in the metadata is set to TRUE if it is found in the consensus
#' It is suggested to run this function _before_ filling the gaps with the function fillGaps()
#'
#' @param segments (cnvstrml)  the segments that will be further annotated
#' @param consensus (cnvstrml) with the consensus
#' @export
#' @return the segments with updated ranges and metadata to reflect if a segment is in consensus
applyConsensus=function(segments, consensus){
  # map them back to each element of the segmentList to get the
  S4Vectors::mcols(segments@segments)['consensus']=NA
  ov=GenomicRanges::findOverlaps( query=segments@segments, subject=consensus@segments)
  #mcols(segments)[ queryHits(ov), 'consensus']=mcols(consensus)[subjectHits(ov),'majority']
  common_segments=GenomicRanges::pintersect( segments@segments[ S4Vectors::queryHits(ov)], consensus@segments[ S4Vectors::subjectHits(ov)])
  ov3=GenomicRanges::findOverlaps( query=common_segments, subject=consensus@segments)
  S4Vectors::mcols(common_segments)[ S4Vectors::queryHits(ov3), 'consensus']=S4Vectors::mcols(consensus@segments)[S4Vectors::subjectHits(ov3),'majority']


  unique_segments=GenomicRanges::setdiff( segments@segments, consensus@segments)
  # add the metadata to the unique_segments
  ov2=findOverlaps( query=unique_segments, subject=segments@segments)
  S4Vectors::mcols( unique_segments )=mcols(segments@segments)[S4Vectors::subjectHits(ov2),]
  segs=c( common_segments, unique_segments) %>% GenomeInfoDb::sortSeqlevels() %>% sort()

  segments@segments=segs
  return(segments)
}


#' get the consensus sequence from a list of GRanges objects that represent the CNV segments
#'
#' The consensus is a new Granges that contains only the shared segments
#' The metadata of the consensus have the cnv.annotations for each of the elements of the segments
#'
#' @param segmentsList (list of GRanges) a list of cnvstrml objects
#' @param method (character) can be either 'qualitative' of 'quantitative'. If it is set to
#' 'qualitative' it returns an object with qualitative information about the segments (e.g NEUT, AMP etc)
#' If it is set to 'quantitative' it retuns an object with the actual number of copy numbers.
#' @returns a Granges with all the segments that are shared across all samples.
#'
findConsensusGR=function( segmentsList, method=c("qualitative","quantitative")){
  # find the regions that are shared across all predictions
  method=match.arg(method)
  # get the shared ranges
  consensus=segmentsList %>%  Reduce( function(x,y){
    ov=GenomicRanges::findOverlaps( query=x, subject=y)
    GenomicRanges::pintersect( x[ S4Vectors::queryHits(ov)], y[ S4Vectors::subjectHits(ov)])

  },.)
  mcols(consensus)=NULL

  # add qualitative information to the consensus
  if(method =='qualitative'){
    metadata=segmentsList %>% lapply( function(X){
      d=rep(NA, length(consensus))
      ov=GenomicRanges::findOverlaps( query=consensus, subject=X)
      d[queryHits(ov)]=X[S4Vectors::subjectHits(ov)]$cnv.annotation
      d
    } ) %>% dplyr::bind_cols()

    majority=metadata %>% apply(. , 1, function(x){ aa=sort(table(x),decreasing = TRUE); names(aa)[1]})
  }

  # add quantitative information to the consensus
  if(method =='quantitative'){
    metadata=segmentsList %>% lapply( function(X){
      d=rep(NA, length(consensus))
      ov=GenomicRanges::findOverlaps( query=consensus, subject=X)
      d[queryHits(ov)]=X[S4Vectors::subjectHits(ov)]$total.allele.copies
      d
    } ) %>% dplyr::bind_cols()
    majority=metadata %>% apply(. , 1, function(x){ median(x,na.rm=TRUE)})
  }

  if( !is.null(names(segmentsList))){
    colnames(metadata)=names( segmentsList)
  }
  S4Vectors::mcols(consensus)=metadata

  S4Vectors::mcols(consensus)$majority=majority
  consensus
}



#' Load a ControlFreec _CNV output from a tab delimited file
#' this function can also use a _CNV file with pvalues added by the accompanying script in the controlFreec package.
#' it expects to find a tab delimited file
#' The columns are:
#' * chromosome
#' * start position
#' * end position
#' * copy.number
#' * status
#' * genotype
#' * uncertainty of the assignment
#' * WilcoxonRankSumTestPvalue
#' * KolmogorovSmirnovPvalue
#' Note: at least for v 11.5 it is not possible to understand how subclones are coded in the output files.
#' As a result we only parse the main clone information (i.e no subclones). A related ticket has been
#' opened on GitHub to help understand the C-FREEC outputs (https://github.com/BoevaLab/FREEC/issues/84)
#' @param cfreec.fn input filename with CNV segments,
#' @param filterLowUncerntainty if set to TRUE segments with uncertainty equal to -1 will be discarded
#' @param pvalThreshold  threshold of pvalues to filter if the pvalues are included in the file
#' @return A GRanges object with the segments predicted by ControlFreec
#' @export
parseCfreec_file=function(cfreec.fn){
  controlfreec=read.table(cfreec.fn,header=FALSE, sep="\t")
  if( controlfreec[1,1]=="chr"){
    colnames(controlfreec)=controlfreec[1,]
    controlfreec=controlfreec[-1,]
  }else{
    colnames(controlfreec)=c( "chr","start" ,"end","copy.number","status","genotype", "uncertainty" )
  }
  return(parseCfreec( controlfreec ))
}

#' Generate a cnv Granges from a dataframe created by ControlFreec (file ending in .CNVs or contains pvalues.)
#' Control-FREEC is a tool for detection of copy-number changes and allelic imbalances (including LOH)
#' using deep-sequencing data. It automatically computes, normalizes, segments copy number and beta allele
#' frequency (BAF) profiles, then calls copy number alterations and LOH. The control (matched normal) sample
#' is optional for whole genome sequencing data but mandatory for whole exome or targeted sequencing data.
#' Control-FREEC up to now (Jan 2021) is the standard tool used in the hCELG and hBMS pipelines for data
#' analysis of WES and WGS datasets.
#' When the genotype is - it means that the progrm did not find any SNPs to decide the correct genotype, and
#' the estimated copy number depends on the logRR
#' The output includes
#' * chromosome
#' * start position
#' * end position
#' * copy.number
#' * status
#' * genotype
#' * uncertainty of the assignment
#' * WilcoxonRankSumTestPvalue
#' * KolmogorovSmirnovPvalue
#' Note: at least for v 11.5 it is not possible to understand how subclones are coded in the output files.
#' @param controlfreec.df input filename with CNV segments,
#' @param filterLowUncerntainty if set to TRUE segments with uncertainty equal to -1 will be discarded
#' @param pvalThreshold  threshold of pvalues to filter if the pvalues are included in the file
#' @return A GRanges object with the segments predicted by ControlFreec
#' @export
parseCfreec=function(controlfreec.df,
                     filterLowUncertainty=TRUE,
                     pvalThreshold=0.01){

  controlfreec.gr=GenomicRanges::makeGRangesFromDataFrame( controlfreec.df,
                                                           keep.extra.columns = TRUE,
                                                           ignore.strand = TRUE)%>%
    GenomeInfoDb::sortSeqlevels( ) %>%
    BiocGenerics::sort()

  GenomeInfoDb::seqlevelsStyle( controlfreec.gr)='UCSC'

  mcols_new=S4Vectors::mcols( controlfreec.gr ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate( total.allele.copies=as.numeric(copy.number),
                   major.allele.copies=stringr::str_count( genotype, 'A'),
                   minor.allele.copies=stringr::str_count( genotype, 'B'),
                   uncertainty=as.numeric(uncertainty)
    )
  flags=apply(mcols_new, 1, function(K){ cnvFlag(total.cn = K['total.allele.copies'], minor.cn = K['minor.allele.copies'] )}) # pass the total.cn to make sure that cases with uncertainty=-1 are parsed properly
  mcols_new = mcols_new %>%
    dplyr::mutate( cnv.flag=flags)
  S4Vectors::mcols( controlfreec.gr)=mcols_new
  if( filterLowUncertainty==TRUE){
    keep=which(controlfreec.gr$uncertainty >-1)
    controlfreec.gr=controlfreec.gr[keep]
  }
  if( "WilcoxonRankSumTestPvalue" %in% colnames( S4Vectors::mcols(controlfreec.gr))){
    keep=which( controlfreec.gr$WilcoxonRankSumTestPvalue< pvalThreshold &&	controlfreec.gr$KolmogorovSmirnovPvalue<pvalThreshold)
    controlfreec.gr=controlfreec.gr[keep]
  }

  return(controlfreec.gr)
}

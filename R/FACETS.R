

#' create a cnvstrml object from a dataframe with facets data
#'
#' @param df (data.frame)  dataframe produced by FACETS
#' @param purity (numeric)  the purity produces by FACETS (a value between 0 - 1)
#' @param ploidy (numeric) the ploidy produced by FACETS
#' @param genome (character) the genome version assumed for this dataset
#' @param sample (character)  the name of the sample
#' @return a cnvstrml object
#'
parseFACETS_df=function(
    df,
    purity=NA_real_,
    ploidy=NA_real_ ,
    sample=NA_character_,
    genome='hg38'
){

  # convert the input to a dataframe
  if(! is.data.frame(df)){
    df=as.data.frame(df)
  }
  # FACETS tend to call Chr X as chr23
  df[ which(df$chrom==23),"chrom"]="X"

  # convert to GRanges
  # and normalize the chromosome name to UCSC
  facets.gr=GenomicRanges::makeGRangesFromDataFrame( df,
                                                     keep.extra.columns = TRUE,
                                                     ignore.strand = TRUE) %>%
    GenomeInfoDb::sortSeqlevels( ) %>%
    BiocGenerics::sort()
  GenomeInfoDb::seqlevelsStyle(facets.gr)='UCSC'


  # Add the flags

  mcols_new=S4Vectors::mcols( facets.gr ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate( total.allele.copies=tcn.em,
                   major.allele.copies=ifelse( is.na(lcn.em), tcn.em, tcn.em-lcn.em),
                   minor.allele.copies=lcn.em,
                   cnv.ccf=cf.em)

  S4Vectors::mcols( facets.gr )=mcols_new

  obj=cnvstrml(
    segments=facets.gr,
    purity=as.numeric(purity),
    ploidy=as.numeric(ploidy),
    method='FACETS',
    sample=sample,
    genome=genome
  )
  obj
}



#' parse a VCF file in BND format with FACETS data
#'
#' Such a file is created from the cnv_facets tool
#'
#' @param vcf (character) the input VCF file
#' @param genome (character)  the genome version assumed for this dataset
#' @param sample (character)  the name of the sample
#' @return a cnvstrml object
#'
parseFACETS_vcf=function(
  vcf,
  sample=NA_character_,
  genome='hg38'
){

  ff.vcf=VariantAnnotation::readVcf( vcf , genome=genome)
  # convert the rowRanges to a dataframe

  ff.gr=SummarizedExperiment::rowRanges(ff.vcf)
  keep=which(ff.gr$FILTER=='PASS')
  ff.vcf=ff.vcf[keep]
  ff.gr=SummarizedExperiment::rowRanges(ff.vcf)

  S4Vectors::mcols(ff.gr)=cbind( S4Vectors::mcols(ff.gr), VariantAnnotation::info( ff.vcf))
  mcols_new=S4Vectors::mcols( ff.gr ) %>%
    tibble::as_tibble() %>%
    tidyr::replace_na(list(LCN_EM=0, CF_EM=0)) %>%
    dplyr::mutate( total.allele.copies=TCN_EM,
                   minor.allele.copies=LCN_EM,
                   major.allele.copies= TCN_EM-LCN_EM,
                   cnv.ccf=CF_EM)

  S4Vectors::mcols( ff.gr )=mcols_new
  GenomicRanges::end(ff.gr)=ff.gr$END
  if(is.na(sample) || length(sample)==0){
    sample=VariantAnnotation::samples(VariantAnnotation::header(ff.vcf))
    #message("Extracting sample name from the header of the VCF file :" , sample)
  }
  if(is.na(sample) || length(sample)==0){
    sample=vcf %>% basename %>% stringi::stri_replace_all_fixed(., c('.vcf.gz','.vcf.bgz','.vcf'), replacement='', vectorize_all = FALSE)
    #message("Extracting sample name from the filename :" , sample)
  }
  obj=cnvstrml(
    segments=ff.gr,
    purity = as.numeric(VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$purity[1,1]),
    ploidy = as.numeric(VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$ploidy[1,1]),
    method='FACETS',
    sample=sample,
    genome=genome
  )
  obj
}

#' parse the output of FACETS
#'
#' this function is a high level parser that can handle
#' either data frames from FACETS or vcf files
#'
#' @param facets.data (character or data.frame) the input data.
#'        if the type of the argument is character it is assumed to be a file.
#'          In this case it expects the extension .vcf or .vcf.gz to parse it as a vcf file.
#'        if the type of the argument is data.frame it is parsed as such.
#' @param genome (character)  the genome version assumed for this dataset
#' @param sample (sample)  the name of the sample. If the input is a vcf file it will use the name of the sample in the file.
#' @return a cnvstrml object
#' @export
parseFACETS=function(
  facets.data,
  sample=NA_character_,
  genome='hg38'
){


  #dispatch to the right function

  if( is.data.frame( facets.data)){
    message("Preparing the segments from a data frame. Purity and Ploidy have to be set independently.")
    obj=parseFACETS_df( facets.data, sample=sample, genome=genome)
  }
  if( is.character(facets.data) && file.exists( facets.data) &&
      (endsWith( facets.data, "vcf") || endsWith(facets.data,"vcf.gz") || endsWith(facets.data, "vcf.bgz"))
    ){
        message("Preparing the segments from the vcf file ", basename( facets.data),". Purity and Ploidy are expected to be found in the header of the file.")
        obj=parseFACETS_vcf( facets.data,sample=sample, genome=genome)
      }

  mcols_new= S4Vectors::mcols(obj@segments) %>% as.data.frame
  flags=apply(mcols_new, 1, function(K){ cnvstrmliner::annotateCNV(major.cn = K['major.allele.copies'], minor.cn = K['minor.allele.copies'] )})

  S4Vectors::mcols(obj@segments) = mcols_new %>%
    dplyr::mutate( cnv.annotation=flags,
                   cnv.loh=ifelse( grepl("LOH",cnv.annotation), TRUE, FALSE),
                   cnv.lrr=CNLR_MEDIAN)

  obj

}

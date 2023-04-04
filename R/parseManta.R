#' Parse the manta ouput (vcf v4.1)
#'
#' @param sv.vcf (character) a file name (in vcf format)
#' @param genome (character) the genome version to use (hg38). This is passed to VariantAnnotation::readVcf
#' @param sample (character) the name of the sample
#' @return a cnvstrml object, Purity and Ploidy are not set
#' @export


parseManta=function(sv.vcf ,sample,genome='hg38'){

  # load the file
  vcf=VariantAnnotation::readVcf( sv.vcf)
  gr=StructuralVariantAnnotation::breakpointRanges(vcf) %>%
       #StructuralVariantAnnotation::breakendRanges(vcf)
    GenomeInfoDb::sortSeqlevels() %>%
    S4Vectors::sort()

  #keep only the svs that pass the filters
  gr=gr[ gr$FILTER == 'PASS']
  obj=cnvstrml(
    segments=gr,
    purity = NA,
    ploidy = NA,
    method='manta',
    sample=sample,
    genome=genome
  )
  obj

}


#" match cnv segments with sv breakpoints
#'
#' In order to match the SV breakpoints to the CNV segments
#' we expect to find a breakpoint in the vicinity of the
#' starting and ending coordinates of the CNV segment
#' Furthermore the two breakpoints need to be partners.
#'
#'
#' @param segments (cnvstrml) the segments predicte to be copy altered
#' @param sv (cnvstrml) the SV ranges from breakpointRanges
#' @param wiggle (numeric) distance between the CNV boundaries and sv breakpoints
#' @export
#' @return an updates GRanges with segments and support from SV
matchCNVSV=function( segs, sv, wiggle=1000){
  cnv=segs@segments
  gr=sv@segments
  end(cnv)=start(cnv)+1
  cnv=flank(cnv, wiggle)
  ov.start=findOverlaps( query=cnv, subject=gr)
  svs.start= data.frame(
          'segment'=names(a@segments)[queryHits(ov.start)],
          'segment.chr'=seqnames(a@segments)[queryHits(ov.start)],
          'segment.start'=start(a@segments)[queryHits(ov.start)],
          'segment.end'=end(a@segments)[queryHits(ov.start)],
          'segment.width'=width(a@segments)[queryHits(ov.start)],
          'segment.type'=S4Vectors::mcols(a@segments)[ queryHits(ov.start),'cnv.annotation'],
          'SV5'=names(gr)[subjectHits(ov.start)],
          'SV5.chr'=seqnames(gr[subjectHits(ov.start)]),
          'SV5.start'=start(gr[subjectHits(ov.start)]),
          'SV5.end'=end(gr[subjectHits(ov.start)]),
          'SV5.len'=mcols(gr)[subjectHits(ov.start),'svLen'],
          'SV5.type'=S4Vectors::mcols(gr)[subjectHits(ov.start), 'svtype'],
          'SV5.partner'=gr[subjectHits(ov.start)]$partner) %>%
    dplyr::mutate( SV5.partner.chr=seqnames(gr[ SV5.partner]) %>% as.vector(),
                   SV5.partner.start=start(gr[ SV5.partner]) ,
                   SV5.partner.end=end(gr[ SV5.partner])
    )

  # svs.start=mergeByOverlaps( query=cnv, subject=gr)
  # svs.start=svs.start[,c('cnv','cnv.annotation', 'cnv.ccf','gr', 'sourceId','partner','svtype' )]
  cnv=segs@segments
  start(cnv)=end(cnv)-1
  cnv=flank(cnv, wiggle)
  ov.end=findOverlaps( query=cnv, subject=gr)
  svs.end = data.frame(
    'segment'=names(a@segments)[queryHits(ov.end)],
    'segment.chr'=seqnames(a@segments)[queryHits(ov.end)],
    'segment.start'=start(a@segments)[queryHits(ov.end)],
    'segment.end'=end(a@segments)[queryHits(ov.end)],
    'segment.width'=width(a@segments)[queryHits(ov.end)],
    'segment.type'=S4Vectors::mcols(a@segments)[ queryHits(ov.end),'cnv.annotation'],
    'SV3'=names(gr)[subjectHits(ov.end)],
    'SV3.chr'=seqnames(gr[subjectHits(ov.end)]),
    'SV3.start'=start(gr[subjectHits(ov.end)]),
    'SV3.end'=end(gr[subjectHits(ov.end)]),
    'SV3.len'=mcols(gr)[subjectHits(ov.end),'svLen'],
    'SV3.type'=S4Vectors::mcols(gr)[subjectHits(ov.end), 'svtype'],
    'SV3.partner'=gr[subjectHits(ov.end)]$partner)%>%
    dplyr::mutate( SV3.partner.chr=seqnames(gr[ SV3.partner]) %>% as.vector(),
                   SV3.partner.start=start(gr[ SV3.partner]) ,
                   SV3.partner.end=end(gr[ SV3.partner])
    )
  #
  # svs.end=mergeByOverlaps( query=cnv, subject=gr)
  # svs.end=svs.end[,c('cnv','cnv.annotation', 'cnv.ccf','gr', 'sourceId','partner','svtype' )]
  #segments where the start and end have a hit to SV and the other end of the segment to the parnter



  rr=dplyr::inner_join( svs.start, svs.end, by=dplyr::join_by( segment==segment)  , multiple='all', suffix=c(".start",".end")) %>%
    dplyr::select( segment,
                   segment.chr=segment.chr.start,
                   segment.start=segment.start.start,
                   segment.end=segment.end.start,
                   segment.width=segment.width.start,
                   segment.type=segment.type.start,
                   SV5,
                   SV5.chr,
                   SV5.start,
                   SV5.end,
                   SV5.type,
                   SV5.partner,
                   SV5.partner.chr,
                   SV5.partner.start,
                   SV5.partner.end,
                   SV5.len,
                   SV3,
                   SV3.chr,
                   SV3.start,
                   SV3.end,
                   SV3.type,
                   SV3.partner,
                   SV3.partner.chr,
                   SV3.partner.start,
                   SV3.partner.end,
                   SV3.len) %>%
    dplyr::mutate( support=1,
                   support=ifelse( SV3.chr==SV5.chr, support+1, support),
                   support=ifelse( SV3.partner==SV5, support+2, support),
                   support=ifelse( segment.type =='AMP' & (SV5.type == 'DEL'|SV3.type=='DEL'), 0 , support),
                   support=ifelse( segment.type =='HOMDEL' & (SV5.type == 'DUP'|SV3.type=='DUP'), 0 , support),
                   support=ifelse( segment.type =='HETDEL' & (SV5.type == 'DUP'|SV3.type=='DUP'), 0 , support)
                    )%>%
    dplyr::filter(support>0)
  #dplyr::inner_join( svs.start, svs.end, by=dplyr::join_by( segment.start==segment.end , SV.start==SV.end.partner, SV.start.partner==SV.end)  , multiple='all')

  S4Vectors::mcols(segs@segments)$sv.support=0

  return(segs)
}



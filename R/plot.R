#' plot the chromosomes with the segments with copy number
#'
#' @param segments (cnvstrml) the object to plot
#' @param max.copies (integer) the maximum number of copies to plot. Any segment with more copies will be downsized to this number. The Y scale of the plot ranges between 0 and max.copies.
#' @export

plotSegments=function( segments, max.copies=6 , main=NULL){
  aa=segments@segments
  if(segments@segtype == 'consensus'){
    message("This function can plot only segments for a single sample")
    return(NULL)
  }
  if( is.null(main)){
    main=segments@sample
  }
  subtitle=paste( segments@method, collapse=",")

  aa %>%
    as.data.frame(row.names = NULL) %>%
    dplyr::mutate( total.allele.copies=ifelse(total.allele.copies > 6, 6,total.allele.copies)) %>%
    ggplot( aes(y=total.allele.copies, x=start)) +
    geom_segment(aes(yend=total.allele.copies, xend=end, color=cnv.annotation), size=4) +
    ylim(0,6) +
    facet_grid(~ seqnames, space="free_x" ,scales='free_x') +
    labs( title=main, subtitle=subtitle ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x=element_blank()) +
    scale_color_manual( values=c('HOMDEL'='blue','HETDEL'='blue','NEUT'='green','NEUT-LOH'='cyan','AMP-LOH'='orange','AMP'='red'))

}



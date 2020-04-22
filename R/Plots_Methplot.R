#' @title Methylation Plot
#'
#' @description This function plots the methylation of preprocessed benchmark and optionally submitted user data
#' for each region in a scatter plot.
#' @usage epiGB.MethPlot(epiGB.data = epiGB.data, user_input = NULL,
#'                       pdf = T, pdf_path = NULL, filename = "epiGB.methplot")
#' @param pdf Optional pdf output, defaults to TRUE.
#' @param pdf_path Optional input of a path where the pdf is saved. Defaults to NULL, the pdf is then saved
#' in the current folder.
#' @param filename Optional change of the core filename, default name is "epiGB.methplot".
#' @return This function generates scatter plots visualizing the methylation levels of all pipelines
#' against all samples, per region. The consensus corridor for each region and sample and its mean is given as an error bar.
#' @import ggplot2
#' @import ggforce
#' @import scico
#' @examples Plotting (without pdf output) the implemented benchmarking data. Note that it is not necessary to supply 'epiGB.data'
#' since it is implemented within the package.
#' epiGB.MethPlot(pdf = F)
#' @export
#'
epiGB.MethPlot <- function(epiGB.data = NULL, user_input = NULL,
                           pdf = T, pdf_path = NULL, filename = "epiGB.methplot", use.cols=F){


  #result of the Preprocessing function is a dataframe called "epiGB.data",
  #including the user data already in the dataframe if entered in the function call
  if (use.cols){
    colors <- ggsci::pal_rickandmorty(palette = "schwifty")(length(unique(epiGB.data$pipeline)))

  #col_list <- c(
    # "methCT-OTP" = colors[1],
    # "methCT-DKFZ_B060" = colors[2],
    # "methCT-BioMedX" = colors[3],
    # "gemBS-IHEC" = colors[4],
    # "bistro" = colors[5],
    # "GSNAP-UdS.NT" = colors[6],
    # "GSNAP-UdS.TT" = colors[7],
    # "LaJolla-WGBS" = colors[8],
    # "DKFZ_B370" = colors[9],
    # "BAT" = colors[10],
    # "gemBS-2" = colors[11])
    col_list <- colors
    names(col_list) <- unique(epiGB.data$pipeline)

  #col_list <- col_list[as.character(unique(epiGB.data$pipeline))]
  }
  methplots <- list()
  epiGB.data$sample <- factor(epiGB.data$sample, levels=c("5T", "6T", "5N", "6N"))


  for (i in c(1:3)){
    methplots[[i]] <- ggplot2::ggplot()+
      geom_point(data=epiGB.data, aes(x = sample, y = Methylation,
                                      size = Coverage*1.2, group = pipeline), colour = "grey10",
                 position=position_dodge(width=0.5))+
      geom_point(data=epiGB.data, aes(x = sample, y = Methylation,
                                      size = Coverage, colour = pipeline),
                 alpha = 0.9, position = position_dodge(width=0.5))+
      geom_errorbar(data=epiGB.data,
                    mapping=aes(x=epiGB.data$sample,
                                ymin=as.numeric(epiGB.data$lower),
                                ymax=as.numeric(epiGB.data$upper)),
                    width=0.2, size=1, color="black")+
      geom_point(data=epiGB.data, aes(x = sample, y = mean, size = 0.5))+
      facet_wrap_paginate(~locus_identifier, nrow = 4, ncol = 4, page = i)+
      labs(title = "DNA Methylation by Region and Sample")+
      labs(colour = "Pipelines")+
      labs(size = "Coverage")+
      guides(colour = guide_legend(override.aes = list(alpha = 1)),
             size = guide_legend(override.aes = list(alpha = 0.1)))+
      ylab("methylation")+
      xlab("samples")+
      scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1))+
      theme_bw()+#this gets rid of the grey background
      coord_flip()+#flip the axes
      theme(strip.background =element_rect(fill="white"))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    if (use.cols){
      methplots[[i]] <-  methplots[[i]] + scale_color_manual(values=col_list, breaks=names(col_list))} else{
        methplots[[i]] <-  methplots[[i]] + scale_color_brewer(palette="Spectral")
      }
    if (pdf==T){
      pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename,".pdf"), width = 13, height = 11)
      invisible(lapply(methplots, print))
      dev.off()
    }
  }
  return(methplots)
}

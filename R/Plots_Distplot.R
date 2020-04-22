#' @title Distance Plots
#'
#' @description This function plots the relative and/or absolute distance of preprocessed benchmark and
#' optionally submitted user data from the consensus corridor either by sample or by region.
#' Plotting (without pdf output) the implemented benchmarking data. Note that it is not necessary to supply 'epiGB.data'
#' since it is implemented within the package.
#' The function allows for the user to choose all or some of the plots to be displayed (all function calls default to TRUE).
#' Generating plots for absolute and relative distance solely by sample:
#' epiGB.DistPlot(by_region = F, sum = F, abs = F, rel = F, pdf = F)
#' @usage epiGB.DistPlot(epiGB.data = epiGB.data, user_input = NULL, by_sample = T, by_region = T, sum = T, abs = T, rel = T,
#'                       pdf=T, pdf_path=NULL, filename="epiGB.distplot" )
#' @param user_input Optional user data input, defaults to NULL.
#' @param by_sample Plots the relative and/or absolute distance by sample as violin plot, defaults to TRUE.
#' @param by_region Plots the relative and/or absolute distance by region as scatter plot, defaults to TRUE.
#' @param sum Plots the average distance over all samples and regions as scatter plot.
#' @param rel Plotts the relative distance for all distance plots, defaults to TRUE.
#' @param abs Plotts the absolute distance for all distance plots, defaults to TRUE.
#' @param pdf Optional pdf output, defaults to TRUE.
#' @param pdf_path Optional input of a path where the pdf is saved. Defaults to NULL, the pdf is then saved
#' in the current folder.
#' @param filename Optional change of the core filename, default name is "epiGB.distplot".
#' @return this function generates a total of six plots (per data set) serving to visualize the distance of methylation measurements
#' from the consensus corridor. Either relative or absolute distance of each pipeline by region as scatter plot as well
#' as by sample in a violin plot is displayed. Additionally, the mean of either absolute or relative distance of each
#' pipeline over all regions is plotted by sample as scatter plot to summarize the accuracy of each pipeline.
#' @import ggplot2
#' @import dplyr
#' @export


epiGB.DistPlot <- function(epiGB.data = NULL, user_input = NULL, by_sample = T, by_region = T, sum = T, abs = T, rel = T, summarized=T,
                           pdf=T, pdf_path=NULL, filename="epiGB.distplot", use.cols=F, height=13, width=11){
  #if (!(is.null(user_input))){
  # DataPreprocessing(user_data = user_input)
  #}else{
  #DataPreprocessing()
  #}
  epiGB.data$sample <- factor(epiGB.data$sample, levels=c("5T", "6T", "5N", "6N"))

  if (use.cols){
    colors <- ggsci::pal_rickandmorty(palette = "schwifty")(length(unique(epiGB.data$pipeline)))

    col_list <- colors
    names(col_list) <- unique(epiGB.data$pipeline)

    col_list2 <- c(
      "5"= colors[3],
      "6" = colors[5]
    )
  } else {
    col_list <- NULL
  }


  if (by_sample == T){
    if (abs == T){
      subset_data_absolute <- list()
      for(i in epiGB.data$sample){
        subset_data_absolute[[i]] <- subset(epiGB.data, sample == i,
                                            select=c(sample, pipeline, abs_dist))
      }

      abs_plotlist <- list()
      for (i in unique(epiGB.data$sample)){
        abs_plotlist[[i]] <- summaryPlot_1(data = subset_data_absolute[[i]],y = "abs_dist",
                                           title = paste0("Absolute Distance from the Consensus Corridor over all Regions for sample ", i),
                                           ylab = "absolute distance from consensus corridor",
                                           breaks = seq(0, 0.5, by = 0.1), limits = c(0,0.5), col_list=col_list)
      }
      if (pdf==T){
        pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_by_sample_absolute", ".pdf"),
            width = width, height = height)
        invisible(lapply(abs_plotlist, print))
        dev.off()
      } else {
        lapply(abs_plotlist, print)

      }
    }

    if (rel == T){
      subset_data_relative <- list()
      for(i in epiGB.data$sample){
        subset_data_relative[[i]] <- subset(epiGB.data, sample == i,
                                            select=c(sample, pipeline, rel_dist))
      }

      rel_plotlist <- list()


      for (i in unique(epiGB.data$sample)){
        rel_plotlist[[i]] <- summaryPlot_1(data = subset_data_relative[[i]],y = "rel_dist",
                                           title = paste0("Relative Distance from the Consensus Corridor over all Regions for Sample", i),
                                           ylab = "relative distance from consensus corridor",
                                           breaks = seq(-0.5, 0.5, by = 0.25), limits = c(-0.5,0.5), col_list=col_list)
      }

      if (pdf==T){
        pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_by_sample_relative", ".pdf"),
            width = width, height = height)
        invisible(lapply(rel_plotlist, print))
        dev.off()
      } else {
        lapply(rel_plotlist, print)
      }
    }
  }





  if (by_region == T){
    if (abs == T){
      abs_dist_plots <- list()


#      for (i in c(1:3)){
        abs_dist_plots <- summaryPlot_2(data = epiGB.data, y = "abs_dist",
                                             title = "Absolute Distance by Regions and Samples",
                                             breaks = seq(0, 0.6, by = 0.2), limits = c(0, 0.6), col_list=col_list)

        if (pdf==T){
          pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_by_region_absolute",
                     ".pdf"), width = width, height = height)
          invisible(lapply(abs_dist_plots, print))
          dev.off()
        } else {
          lapply(abs_dist_plots, print)
        }
#      }
    }

    if (rel == T){

      rel_dist_plots <- list()

      #for (i in c(1:3)){
        rel_dist_plots <- summaryPlot_2(data = epiGB.data, y = "rel_dist",
                                             title = "Relative Distance by Regions and Samples",
                                             breaks = seq(-0.5, 0.5, by = 0.25), limits = c(-0.5, 0.5), col_list=col_list)

        if (pdf==T){
          pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_by_region_relative",
                     ".pdf"), width = width, height = height)
          invisible(lapply(rel_dist_plots, print))
          dev.off()
        } else {
          lapply(rel_dist_plots, print)
        }
      #}
    }
  }

  if (sum == T){
    if(abs == T){
        #mean_by_pipeline <- NULL
      #browser()
      mean_by_pipeline <- epiGB.data %>%
        dplyr::group_by(pipeline, sample) %>%
        dplyr::summarize(av_dist = mean(abs_dist, na.rm=T), av_coverage = mean(Coverage, na.rm=T))
      mean_by_pipeline$patient <- substr(mean_by_pipeline$sample, 1, 1)
      mean_by_pipeline$type <- substr(mean_by_pipeline$sample, 2, 2)

      abs_dist <- summaryPlot_3(data=mean_by_pipeline, x = "pipeline", y = "av_dist", size = "av_coverage", group = "patient", types="type",
                                ylab="absolute distance", title="Avarage Absolute Distance from the Consensus Corridor", col_list=col_list2)

      if (pdf==T){
        print(abs_dist)
        ggsave(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_sum_absolute",
                      ".pdf"), abs_dist, width = width, height = height)
        dev.off()
      } else {
        print(abs_dist)
      }
    }

    if (rel == T){

      mean_by_pipeline <- epiGB.data %>%
        dplyr::group_by(pipeline, sample) %>%
        dplyr::summarize(av_dist = mean(rel_dist, na.rm=T), av_coverage = mean(Coverage, na.rm=T))
      mean_by_pipeline$patient <- substr(mean_by_pipeline$sample, 1, 1)
      mean_by_pipeline$type <- substr(mean_by_pipeline$sample, 2, 2)


      rel_dist <- summaryPlot_3(data=mean_by_pipeline, x = "pipeline", y = "av_dist", size = "av_coverage", group = "patient", types="type",
                                ylab="directional distance", title="Avarage Relative Distance from the Consensus Corridor", col_list=col_list2)

      if (pdf==T){
        print(rel_dist)
          ggsave(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_sum_relative",
                      ".pdf"), rel_dist, width = width, height = height)
        dev.off()
      } else {
        print(rel_dist)
      }
    }
  }

  if (summarized == T){
    if(abs==T){
      #mean_by_pipeline <- NULL
      #browser()
      sum_by_pipeline <- epiGB.data %>%
        dplyr::group_by(pipeline, sample) %>%
        dplyr::summarize(sum_dist = sum(abs_dist, na.rm=T), av_coverage = mean(Coverage, na.rm=T))
      sum_by_pipeline$patient <- substr(sum_by_pipeline$sample, 1, 1)
      sum_by_pipeline$type <- substr(sum_by_pipeline$sample, 2, 2)

      abs_dist <- summaryPlot_4(data=sum_by_pipeline, x = "pipeline", y = "sum_dist", size = "av_coverage", group = "patient", types="type",
                                ylab="absolute distance", title="Summarized Absolute Distance from the Consensus Corridor", col_list=col_list2)

      if (pdf==T){
        print(abs_dist)
        ggsave(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_summarized_absolute",
                      ".pdf"), abs_dist, width = width, height = height)
        dev.off()
      } else {
        print(abs_dist)
      }
    }
    # if (rel == T){
    #   #mean_by_pipeline <- NULL
    #   #browser()
    #   sum_by_pipeline <- epiGB.data %>%
    #     dplyr::group_by(pipeline, sample) %>%
    #     dplyr::summarize(sum_dist = sum(rel_dist, na.rm=T), av_coverage = mean(Coverage, na.rm=T))
    #   sum_by_pipeline$patient <- substr(sum_by_pipeline$sample, 1, 1)
    #   sum_by_pipeline$type <- substr(sum_by_pipeline$sample, 2, 2)
    #
    #   rel_dist <- summaryPlot_4(data=sum_by_pipeline, x = "pipeline", y = "sum_dist", size = "av_coverage", group = "patient", types="type",
    #                             ylab="relative distance", title="Summarized Relative Distance from the Consensus Corridor", col_list=col_list2)
    #
    #   if (pdf==T){
    #     print(rel_dist)
    #     ggsave(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_summarized_relative",
    #                   ".pdf"), rel_dist, width = width, height = height)
    #     dev.off()
    #   } else {
    #     print(rel_dist)
    #   }
    #
    # }
  }
}

#' @import ggplot2
#' @import ggforce
summaryPlot_1 <- function(data,y,title, ylab = "", breaks, limits, col_list){
  g <- ggplot2::ggplot(data = data, aes(x = factor(pipeline,levels = rev(levels(factor(pipeline)))),
                                        y = get(y), fill = pipeline))+
    geom_violin(alpha = 0.7, scale = "width")+
    labs(title = title)+
    stat_summary(fun.y=median, geom="point", size=2)+
    ylab(ylab)+
    xlab("")+ # maybe add "pipelines"?
    guides(fill = FALSE)+
    scale_y_continuous(breaks = breaks, limits = limits)+
    theme_bw()+#this gets rid of the grey background
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), axis.line = element_line())
  if (!is.null(col_list)){
    g <-g +  scale_fill_manual(values=col_list, breaks=names(col_list))} else{
      g <-g +  scale_fill_brewer(palette="Spectral")
    }
  return(g)
}

#' @import ggplot2
#' @import ggforce
summaryPlot_2 <- function(data, y, title = "", breaks, limits, col_list){
  g <- list()
  pages <- ceiling(length(unique(data$locus_identifier))/16)
  for (i in c(1:pages)){
  g[[i]] <- ggplot2::ggplot()+
    geom_point(data=data, aes(x = sample, y = get(y),
                              size = Coverage*1.2, group = pipeline),
               colour = "grey10", position=position_dodge(width=0.5))+
    geom_point(data=data, aes(x = sample, y = get(y),
                              size = Coverage, colour = pipeline),
               alpha = 0.9, position=position_dodge(width=0.5))+
    facet_wrap_paginate(~locus_identifier, nrow = 4, ncol = 4, page = i)+
    labs(title = title)+
    labs(colour = "Pipelines")+
    labs(size = "Coverage")+
    guides(colour = guide_legend(override.aes = list(alpha = 1)),
           size = guide_legend(override.aes = list(alpha = 0.1)))+
    ylab("Methylation")+ # lowercase?
    xlab("samples")+
    scale_y_continuous(breaks = breaks, limits = limits)+
    theme_bw()+#this gets rid of the grey background
    coord_flip()+#flip the axes
    theme(strip.background =element_rect(fill="white"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  if (!is.null(col_list)){
    g[[i]] <- g[[i]] +  scale_color_manual(values=col_list, breaks=names(col_list))} else{
      g[[i]] <- g[[i]] +  scale_color_brewer(palette="Spectral")
    }
  }
  return(g)
}

#' @import ggplot2
#' @import ggforce
summaryPlot_3 <- function(data, x, y, size, group, types, ylab="", title="", col_list){

  g <- ggplot2::ggplot() + geom_point(data=data, aes(x = get(x), y = get(y),
                                            size = get(size)*1.2, group = get(group), shape=get(types)),
                             colour = "grey40", alpha = 0.7)+
    #,
                             #position=position_dodge(width=0.5))+
    geom_point(data=data, aes(x = get(x), y = get(y),
                              colour = get(group), size = get(size), shape=get(types)),
               alpha = 0.7)+
               #, position=position_dodge(width=0.5))+

    labs(title = title)+
    labs(colour = "Patient")+
    labs(size = "Coverage")+
    labs(shape="Sample Type")+
    guides(colour = guide_legend(override.aes = list(alpha = 1)),
           size = guide_legend(override.aes = list(alpha = 0.1)))+
    ylab(ylab)+
    xlab("Pipelines")+
    theme_bw()+#this gets rid of the grey background
    coord_flip()+#flip the axes
    theme(strip.background =element_rect(fill="white"))+
    theme(panel.grid.minor = element_blank())
  if (!is.null(col_list)){
    g <- g +  scale_color_manual(values=col_list, breaks=names(col_list))} else{
      g <- g +  scale_color_brewer(palette="Spectral")
    }

  return(g)
}
summaryPlot_4 <- function(data, x, y, size, group, types, ylab="", title="", col_list){

  g <- ggplot2::ggplot() + geom_point(data=data, aes(x = get(x), y = get(y),
                                                     size = get(size)*1.2, group = get(group), shape=get(types)),
                                      colour = "grey40", alpha = 0.7)+
                                      #,
                                      #position=position_dodge(width=0.5))+
    geom_point(data=data, aes(x = get(x), y = get(y),
                              colour = get(group), size = get(size), shape=type),
               alpha = 0.7)+
               #, position=position_dodge(width=0.5))+
    labs(title = title)+
    labs(colour = "Samples")+
    labs(size = "Coverage")+
    labs(shape = "Sample Type")+
    guides(colour = guide_legend(override.aes = list(alpha = 1)),
           size = guide_legend(override.aes = list(alpha = 0.1)))+
    ylim(0, max(data[,y], na.rm=T)+0.5)+
    ylab(ylab)+
    xlab("Pipelines")+
    theme_bw()+#this gets rid of the grey background
    coord_flip()+#flip the axes
    theme(strip.background =element_rect(fill="white"))+
    theme(panel.grid.minor = element_blank())
  if (!is.null(col_list)){
    g <- g +  scale_color_manual(values=col_list, breaks=names(col_list))} else{
      g <- g +  scale_color_brewer(palette="Spectral")
    }

  return(g)
}


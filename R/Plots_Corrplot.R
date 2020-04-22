#' @title Correlation Plots
#'
#' @description This function plots the distance of preprocessed benchmark and optionally submitted user data
#' from the consensus corridor as a correlation scatter plot with regression line.
#' Plotting (without pdf output) the implemented benchmarking data. Note that it is not necessary to supply 'epiGB.data'
#' since it is implemented within the package.
#' epiGB.CorrPlot(pdf = F)
#' The function allows for the user to choose all or some of the plots to be displayed (all function calls default to TRUE).
#' Plotting only correlation of the relative distance from the consensus corridor for each sample:
#' @usage epiGB.CorrPlot(epiGB.data = NULL, user_input = NULL, pdf=T, pdf_path=NULL,
#'                       rel = T, meth = T, abs = T, filename="epiGB.correlation")
#' @param user_input Optional user data input, defaults to NULL.
#' @param rel Plotting of the relative distance,  defaults to TRUE.
#' @param meth Plotting of the methylation values, defaults to TRUE.
#' @param abs Plotting of the absolute distance,  defaults to TRUE.
#' @param pdf Optional pdf output, defaults to TRUE.
#' @param pdf_path Optional input of a path where the pdf is saved. Defaults to NULL, the pdf is then saved
#' in the current folder.
#' @param filename Optional change of the core filename, default name is "epiGB.correlation".
#' @return This function generates a matrix of pairwise (between the respective pipelines) correlation scatter plots
#' (for absolute and relative distance from a predefined consensus corridor, as well as methylation levels measured
#' by different pipelines) alongside the respective correlation coefficients (Pearson correlation) and density plots
#' for each sample supplied by the benchmarking (and optionally user) data. If pdf output is requested the function
#' returns three pdf files with the filename extensions _absolute, _methylation and _relative, respectively.
#' @import GGally
#' @import ggplot2
#' @import dplyr
#' @export


epiGB.CorrPlot <- function(epiGB.data = NULL, user_input = NULL, rel = T, meth = T, abs = T,
                           pdf=T, pdf_path=NULL, filename="epiGB.correlation"){
  #if (!(is.null(user_input))){
  #DataPreprocessing(user_data = user_input)
  #}else{
  #DataPreprocessing()
  #}

  regression_line <- function(data, mapping, method="lm", ...){
    p <- ggplot2::ggplot(data = data, mapping = mapping) +
      geom_point() +
      geom_smooth(method=method, ...)
    p
  }

  region_names <- unique(epiGB.data$locus_identifier)
  pipeline <- unique(epiGB.data$pipeline)
  sample <- unique(epiGB.data$sample)
  #browser()
  if (rel==T){
    rel_distance_list <- list()


    for (pip in pipeline){
      rel_distance_list[[pip]] <- epiGB.data[grepl(pip, epiGB.data$pipeline), c("locus_identifier","sample", "rel_dist")]
      colnames(rel_distance_list[[pip]]) <- c("locus_identifier","sample", pip)
      if (length(rel_distance_list)==1){
        rel_distance_values <- rel_distance_list[[pip]]
      } else {
        rel_distance_values <- dplyr::left_join(rel_distance_values, rel_distance_list[[pip]], by=c("locus_identifier", "sample"))
      }
    }




    rel_graph_list <- list()

    for (i in sample){
      df_subset <- subset(rel_distance_values, sample == i,
                          select=c(pipeline))
      df_subset <- rel_distance_values[rel_distance_values$sample == i,pipeline]
      rel_graph_list[[i]] <- GGally::ggpairs(df_subset[1:ncol(df_subset)], lower = list(continuous = regression_line),
                                             title = paste0("Correlation of the Relative Distance from the Consensus Corridor in Sample", " ", i))
    }

    if (pdf==T){
      pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_relative", ".pdf"), width = 13, height = 11)
      invisible(lapply(rel_graph_list, print))
      dev.off()
    } else  {
      lapply(rel_graph_list, print)
    }
  }

  if(meth == T){
    methylation_list <- list()

    for (pip in pipeline){
      methylation_list[[pip]] <- epiGB.data[grepl(pip, epiGB.data$pipeline), c("locus_identifier","sample", "Methylation")]
      colnames(methylation_list[[pip]]) <- c("locus_identifier","sample", pip)
      if (length(methylation_list)==1){
        distance_values <- methylation_list[[pip]]
      } else {
        distance_values <- dplyr::left_join(distance_values, methylation_list[[pip]], by=c("locus_identifier", "sample"))
      }
    }




    meth_graph_list <- list()

    for (i in sample){
      df_subset <- subset(distance_values, sample == i,
                          select=c(pipeline))
      df_subset <- distance_values[distance_values$sample == i,pipeline]
      meth_graph_list[[i]] <- GGally::ggpairs(df_subset[1:ncol(df_subset)], lower = list(continuous = regression_line),
                                              title = paste0("Methylation values", " ", i))
    }

    if (pdf==T){
      pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_methylation", ".pdf"), width = 13, height = 11)
      invisible(lapply(meth_graph_list, print))
      dev.off()
    } else {
      lapply(meth_graph_list, print)
    }
  }

  if(abs == T){
    abs_distance_list <- list()


    for (pip in pipeline){
      abs_distance_list[[pip]] <- epiGB.data[grepl(pip, epiGB.data$pipeline), c("locus_identifier","sample", "abs_dist")]
      colnames(abs_distance_list[[pip]]) <- c("locus_identifier","sample", pip)
      if (length(abs_distance_list)==1){
        distance_values <- abs_distance_list[[pip]]
      } else {
        distance_values <- dplyr::left_join(distance_values, abs_distance_list[[pip]], by=c("locus_identifier", "sample"))
      }
    }


    regression_line <- function(data, mapping, method="lm", ...){
      p <- ggplot2::ggplot(data = data, mapping = mapping) +
        geom_point() +
        geom_smooth(method=method, ...)
      p
    }


    abs_graph_list <- list()

    for (i in epiGB.data$sample){
      df_subset <- subset(distance_values, sample == i,
                          select=c(pipeline))
      df_subset <- distance_values[distance_values$sample == i,pipeline]
      abs_graph_list[[i]] <- GGally::ggpairs(df_subset[1:ncol(df_subset)], lower = list(continuous = regression_line),
                                             title = paste0("Correlation of the Absolute Distance from the Consensus Corridor in Sample", " ", i))
    }

    if (pdf==T){
      pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_absolute", ".pdf"), width = 13, height = 11)
      invisible(lapply(abs_graph_list, print))
      dev.off()
    } else {
      lapply(abs_graph_list, print)
    }
  }
}

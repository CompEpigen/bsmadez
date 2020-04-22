#' @title Heatmaps
#'
#' @description This function creates a heatmap based on the methylation or the relative and/or absolute distance
#' of preprocessed benchmark and optionally submitted user data from the consensus corridor.
# Add this to the examples!!
#' Plotting (without pdf output) the implemented benchmarking data. Note that it is not necessary to supply 'epiGB.data'
#' since it is implemented within the package.
#' epiGB.Heatmap(pdf = F)
#' Plotting of benchmarking data exclusively on the methylation heatmap (without pdf output):
#' @usage epiGB.Heatmap(epiGB.data=epiGB.data, user_input = NULL, meth = T, dist = T, abs = T, rel = T,
#'                      pdf=T, pdf_path=NULL, filename="epiGB.heatmap")
#' @param user_input Optional user data input, defaults to NULL.
#' @param meth Plots the correlation of methylation values of all pipelines.
#' @param dist Plots the correlation of relative and/or absolute distance of all pipelines.
#' @param abs Plots the correlation of absolute distance in both the methylation and distance plot,
#' defaults to TRUE.
#' @param rel Plots the correlation of relative distance in both the methylation and distance plot,
#' defaults to TRUE.
#' @param pdf Optional pdf output, defaults to TRUE.
#' @param pdf_path Optional input of a path where the pdf is saved. Defaults to NULL,
#' the pdf is then saved in the current folder.
#' @param filename Optional change of the core filename, default name is "epiGB.heatmap".
#' @return this function returns two heatmaps visualizing methylation levels and their respective absolute
#' and relative distance from the consensus corridor by sample over all preselected regions. The user decides which
#' of the two heatmaps (methylation and/or relative distance) should be plotted (all function calls default to TRUE).
#' @import data.table
#' @import corrplot
#' @import RColorBrewer
#' @export


epiGB.Heatmap <- function(epiGB.data=NULL, meth = T, dist = T, abs = T, rel = T,
                          pdf=T, pdf_path=NULL, filename="epiGB.heatmap"){



  merged_region_values <- as.data.frame(data.table::dcast(data.table::setDT(epiGB.data),
                                                          chr+start+locus_identifier+sample~method+pipeline, value.var= c("Methylation", "Coverage")))

  samples <- unique(epiGB.data$sample)
  pipelines <- unique(epiGB.data$pipeline)
  method <- unique(epiGB.data$method)
  tb_excluded <- colnames(merged_region_values)[!colnames(merged_region_values) %in% c("chr", "start", "locus_identifier", "sample")]
  meth_excluded <- as.character(tb_excluded[-grep("Methylation_", tb_excluded)])
  cov_excluded <- as.character(tb_excluded[-grep("Coverage_", tb_excluded)])

  if (meth == T){
    heatmap_meth_list <- list()
    cor_list <- list()

    for(i in samples){
      heatmap_meth_list[[i]] <- merged_region_values[grepl(pattern = i,
                                                           merged_region_values$sample), cov_excluded]
      heatmap_meth_list[[i]] <- as.matrix(heatmap_meth_list[[i]])
      colnames(heatmap_meth_list[[i]]) <- gsub(paste0("Methylation_", method, "_"), "", colnames(heatmap_meth_list[[i]]))
      class(heatmap_meth_list[[i]]) <- "numeric"
      # try_catch(expr = x,
      #           .e = ~ paste0("There is an error: ", .x),
      #           .w = ~ paste0("This is a warning: ", .x))
      cor_list[[i]] <- matrix(nrow = length(pipelines), ncol = length(pipelines))
      colnames(cor_list[[i]]) <- as.character(pipelines)
      rownames(cor_list[[i]]) <- as.character(pipelines)

    }

    #for the correlation computation, the values must be numeric
    #create an empty matrix for the correlation values

    #loop over the matrix to compute the correlation and put the output in the new matrix

    for (i in samples){
      for (j in pipelines){
        for (m in pipelines){
          cor_list[[i]][j,m] <- cor(heatmap_meth_list[[i]][,j], heatmap_meth_list[[i]][,m], method = c("pearson"), use="pairwise.complete.obs")
        }
      }
    }

    #plot the methylation heatmap

    meth_heatmaps <- list()

    if (pdf==T){
      pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_methylation", ".pdf"),
          width = 13, height = 11)}
    for (i in samples){
      meth_heatmaps[[i]] <- corrplot::corrplot(cor_list[[i]], type = "upper",
                                               title = paste0("Correlation of Methylation in Sample", " ", i),
                                               tl.col="black",
                                               cl.length = 5,
                                               number.digits = 3,
                                               col = colorRampPalette(brewer.pal(11,"Spectral"))(100),
                                               mar = c(2, 0, 4, 0))
      #cl.lim=c(0.971, 1), is.corr=FALSE,
    }
    if (pdf == T){
      dev.off()
    }
  }

  if (dist == T){
    # if(length(unique(pipelines)) < length(pipelines)){
    #   stop("input pipelines overlap with benchmark pipelines")
    # }
    if(abs == T){
      abs_heatmap_dist_list <- list()
      abs_casted_dist_list <- list()
      abs_cor_dist_list <- list()

      for(i in samples){
        abs_heatmap_dist_list[[i]] <- epiGB.data[grepl(pattern = i, epiGB.data$sample),
                                                 c("pipeline", "abs_dist")]
        abs_heatmap_dist_list[[i]] <- within(abs_heatmap_dist_list[[i]], {
          pipeline <- as.character(pipeline)
          ID <- ave(pipeline, pipeline, FUN = seq_along)
        })

        abs_casted_dist_list[[i]] <- reshape2::dcast(abs_heatmap_dist_list[[i]], ID ~ pipeline, value.var = "abs_dist")
        abs_casted_dist_list[[i]] <- as.matrix(abs_casted_dist_list[[i]])
        class(abs_casted_dist_list[[i]]) <- "numeric"
        abs_cor_dist_list[[i]] <- matrix(nrow = length(pipelines), ncol = length(pipelines))
        colnames(abs_cor_dist_list[[i]]) <- as.character(pipelines)
        rownames(abs_cor_dist_list[[i]]) <- as.character(pipelines)

      }



      #loop over the matrix to compute the correlation and put the output in the new matrix

      for (i in samples){
        for (j in pipelines){
          for (m in pipelines){
            abs_cor_dist_list[[i]][j,m] <- cor(abs_casted_dist_list[[i]][,j], abs_casted_dist_list[[i]][,m],
                                               method = c("pearson"), use="pairwise.complete.obs")
          }
        }
      }
      #plot the distance heatmap
      abs_dist_heatmaps <- list()
      if (pdf == T){
        pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_absolute", ".pdf"), width = 13, height = 11)
      }

      for (i in samples){
        abs_dist_heatmaps[[i]] <- corrplot::corrplot(abs_cor_dist_list[[i]], type = "upper",
                                                     title = paste0("Correlation of Absolute Distance from the Consensus Corridor in Sample", " ", i),
                                                     tl.col="black",
                                                     cl.length = 5,
                                                     number.digits = 3,
                                                     col = colorRampPalette(brewer.pal(11,"Spectral"))(100),
                                                     mar = c(2, 0, 4, 0))
        #cl.lim=c(-0.5, 1),
        # mar = c(2, 0, 2, 0))
      }

      if (pdf==T){
        dev.off()
      }
    }

    if(rel==T){

      rel_heatmap_dist_list <- list()
      rel_casted_dist_list <- list()
      rel_cor_dist_list <- list()
      for(i in samples){
        rel_heatmap_dist_list[[i]] <- epiGB.data[grepl(pattern = i, epiGB.data$sample),
                                                 c("pipeline", "rel_dist")]
        rel_heatmap_dist_list[[i]] <- within(rel_heatmap_dist_list[[i]], {
          pipeline <- as.character(pipeline)
          ID <- ave(pipeline, pipeline, FUN=seq_along)
        })

        rel_casted_dist_list[[i]] <- reshape2::dcast(rel_heatmap_dist_list[[i]], ID ~ pipeline, value.var="rel_dist")
        rel_casted_dist_list[[i]] <- as.matrix(rel_casted_dist_list[[i]])
        class(rel_casted_dist_list[[i]]) <- "numeric"
        rel_cor_dist_list[[i]] <- matrix(nrow = length(pipelines), ncol = length(pipelines))
        colnames(rel_cor_dist_list[[i]]) <- as.character(pipelines)
        rownames(rel_cor_dist_list[[i]]) <- as.character(pipelines)
      }
      #loop over the matrix to compute the correlation and put the output in the new matrix
      for (i in samples){
        for (j in pipelines){
          for (m in pipelines){
            rel_cor_dist_list[[i]][j,m] <- cor(rel_casted_dist_list[[i]][,j], rel_casted_dist_list[[i]][,m],
                                               method = c("pearson"), use = "pairwise.complete.obs")
          }
        }
      }
      #plot the distance heatmap
      rel_dist_heatmaps <- list()
      if (pdf == T){
        pdf(paste0(ifelse(is.null(pdf_path), "./", pdf_path), filename, "_relative", ".pdf"), width = 13, height = 11)}
      for (i in samples){
        rel_dist_heatmaps[[i]] <- corrplot::corrplot(rel_cor_dist_list[[i]], type = "upper",
                                                     title = paste0("Correlation of Relative Distance from the Consensus Corridor in Sample", " ", i),
                                                     tl.col="black",
                                                     cl.length = 5,
                                                     number.digits = 3,
                                                     col = colorRampPalette(brewer.pal(11,"Spectral"))(100),
                                                     mar = c(2, 0, 4, 0))
        #cl.lim=c(-0.5, 1),
        #mar = c(2, 0, 2, 0))
      }
      if (pdf == T){
        dev.off()
      }
    }
  }
}

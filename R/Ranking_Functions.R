# ties.method: think of how to deal with similar ranking values

#' @title function to run rankings based on individual samples
#'
#' @description rankings by sample ------
#' @param epiGB.data
#' @import dplyr
#' @export
rank.bySample <- function(epiGB.data = NULL){

  methods <- names(epiGB.data)
  ranking_list <- list()

  for (method in methods){

    ranking_list[[method]] <- epiGB.data[[method]] %>%
      dplyr::group_by(method, sample, pipeline) %>%
      dplyr::summarise(sum_abs_dist = sum(abs_dist, na.rm = T))
    ranking_list[[method]] <- dplyr::mutate(ranking_list[[method]],
                                     ranking = rank(x = sum_abs_dist, ties.method = c("min"))) # change ties.method?

  }
  ranking.df <- do.call(rbind.data.frame, ranking_list)
  # order by ranking
  ranking.df <- arrange(ranking.df, method, sample, ranking)
  return(ranking.df)
}


# Think of whether using sum() or mean() in the summarise function
# ties.method: think of how to deal with similar ranking values

#' @title function to run rankings based on sample type
#'
#' @description rankings by type ------
#' @param epiGB.data
#' @import dplyr
#' @export
rank.byType <- function(ranking_by_sample = NULL, epiGB.data = NULL){

  ranking_by_sample$type <- substr(ranking_by_sample$sample, 2,2)
  byType <- ranking_by_sample
  byType <-  byType %>%
    group_by(method, type, pipeline) %>%
    summarise(ranking = mean(ranking, na.rm = T))
  # change sum() to mean()? pipeline "BioMedX" was not used on sample 5T in method "TWGBS" --> has one ranking less
  byType <- mutate(byType,
                   ranking = rank(x = ranking, ties.method = c("min"))) # !ties.method
  byType <- arrange(byType, method, type, ranking)
  return(byType)
}


# Think of whether using sum() or mean() in the summarise function
# ties.method: think of how to deal with similar ranking values

#' @title function to run rankings based on method
#'
#' @description rankings by method ------
#' @param epiGB.data
#' @import dplyr
#' @export
#'
rank.byMethod <- function(ranking_by_sample = NULL, epiGB.data = NULL){

  byMethod <- ranking_by_sample
  byMethod <- byMethod %>%
    group_by(method, pipeline) %>%
    summarise(ranking = mean(ranking, na.rm = T )) # sample 5N in TWGBS has one pipeline less than the others therefore use mean instead of sum
  byMethod <- mutate(byMethod,
                     ranking = rank(x = ranking, ties.method = c("min"))) # ties.method
  byMethod <- arrange(byMethod, method, ranking)
  return(byMethod)
}

#' @title function to summarise rankings by sample, type and method and give the best 3 pipelines for each category
#' @description summary ranking ------
#' @param ranking_by_sample
#' @param ranking_by_type
#' @param ranking_by_method
#' @param epiGB.data
#' @import dplyr
#' @import plyr
#' @export
#'
rank.summary <- function(ranking_by_sample = NULL, ranking_by_type = NULL, ranking_by_method = NULL, epiGB.data = NULL){

  ### best 3 pipelines by sample:
  sum_1 <- ranking_by_sample %>%
    group_by(method, sample, ranking) %>%
    select(method, sample, pipeline, ranking) %>%
    filter(ranking %in% 1:3) %>%
    mutate(category = paste0("sample", " ", sample)) %>%
    ungroup() %>%
    select(method, category, ranking, pipeline)

  ### best 3 pipelines by type:
  sum_2 <-  ranking_by_type %>%
    group_by(method, type, ranking) %>%
    filter(ranking %in% 1:3) %>%
    mutate(category = paste0("type", " ",type)) %>%
    ungroup() %>%
    select(method, category, ranking, pipeline)

  ### best 3 pipelines by method:
  sum_3 <- ranking_by_method %>%
    group_by(method, ranking) %>%
    filter(ranking %in% 1:3) %>%
    mutate(category = paste0("method", " ", method)) %>%
    ungroup() %>%
    select(method, category, ranking, pipeline)

  total_sum <- rbind(sum_3, sum_2, sum_1) %>% group_by(method)

  methods <- names(epiGB.data)
  method_list <- list()

  for (method in methods){
    total_sum_method <- total_sum[total_sum$method == method,]
    method_list[[method]] <- total_sum_method
  }

  final_sum_list <- list()
  for (method in methods){

    pipeline_1.df <- method_list[[method]][method_list[[method]]$ranking == 1,] %>%
      rename("pipeline_1" = "pipeline") %>%
      ungroup() %>%
      select(method, category, pipeline_1)
    pipeline_1.df <- aggregate(pipeline_1.df["pipeline_1"], by = c(pipeline_1.df["method"], pipeline_1.df["category"]), paste0)

    pipeline_2.df <- method_list[[method]][method_list[[method]]$ranking == 2,] %>%
      rename("pipeline_2" = "pipeline") %>%
      ungroup() %>%
      select(category, pipeline_2)
    pipeline_2.df <- aggregate(pipeline_2.df["pipeline_2"], by=pipeline_2.df["category"], paste0)

    pipeline_3.df <- method_list[[method]][method_list[[method]]$ranking == 3,] %>%
      rename("pipeline_3" = "pipeline") %>%
      ungroup() %>%
      select(category, pipeline_3)
    pipeline_3.df <- aggregate(pipeline_3.df["pipeline_3"], by=pipeline_3.df["category"], paste0)

    final_sum <- left_join(pipeline_1.df, pipeline_2.df, by = "category",  match = "all") # check type and match options
    final_sum <- left_join(final_sum, pipeline_3.df, by = "category",  match = "all") # "

    final_sum_list[[method]] <- final_sum
  }

  final_sum.df <- ldply (final_sum_list, data.frame, .id = NULL)
  # final_sum.df_test <- do.call(rbind, lapply(final_sum_list, as.data.frame)) # better because it doesn't require plyr?
  return(final_sum.df)
}


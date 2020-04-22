# Ranking.Master

#' @title Main function to run all rankings and summarise the results
#'
#' @description All rankings ------
#' @param benchmarking_data A methrix object containing the bechmarking data. If not provided, the package's own data will be used.
#' @param user_data Optional user data input in a methrix object. The sample annotation should have the following columns: sample, method, assembly and pipeline.
#' The samples available are: 5T, 5N, 6T, 6N. The methods can be: WGBS, TWGBS, PBAT or SWIFT. The assembly has to be hg38.
#' @param wdir The output directory. If not provided, the working directory will be used. If non-existing, it will be created.
#' @export


rankings <- function(benchmarking_data=NULL, user_data=NULL, verbose=TRUE, wdir=NULL){

  if (is.null(wdir)){
    if (verbose)
      cat("Using the default working directory as an output directory. /n")
    wdir <- getwd()
  } else if (dir.exists(wdir)){
    if (verbose)
      cat("Using the ", wdir, " as an output directory. /n")
  } else {
    dir.create(file.path(wdir))
    cat("Directory ", wdir, "doesn't exist. Creating it and using it as an output directory. /n")
  }

  if (is.null(benchmarking_data)){
    if (verbose)
      cat("No benchmarking data was provided. Using the package's own.")
    benchmarking_data <- bsmadez::benchmarking_data
  }

# option for user data?
  # data_to_use <- list("benchmarking_data" = benchmarking_data)
  # if (!is.null(user_data)){
  #   data_to_use <- list("benchmarking_data" = benchmarking_data, "user_data" = user_data)
  # }

df <- bsmadez::Preprocess.MethData(benchmarking_data)

# maybe if user_data exists, put benchmarking_data and user_data in a list
# then iterate around the list

ranking_by_sample <- rank.bySample(epiGB.data = df)
ranking_by_type <- rank.byType(ranking_by_sample = ranking_by_sample, epiGB.data = df)
ranking_by_method <- rank.byMethod(ranking_by_sample = ranking_by_sample, epiGB.data = df)
ranking_summary <- rank.summary(ranking_by_sample = ranking_by_sample, ranking_by_type = ranking_by_type,
                                ranking_by_method = ranking_by_method, epiGB.data = df)

ranking_list <- list("ranking_by_sample" = ranking_by_sample, "ranking_by_type" = ranking_by_type,
                     "ranking_by_method" = ranking_by_method, "ranking_summary" = ranking_summary)
return(ranking_list)
}

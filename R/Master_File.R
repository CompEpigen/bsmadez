#' @title Main function to run all plots and analysis
#'
#' @description All plots ------
#' @param benchmarking_data A methrix object containing the bechmarking data. If not provided, the package's own data will be used.
#' @param user_data Optional user data input in a methrix object. The sample annotation should have the following columns: sample, method, assembly and pipeline.
#' The samples available are: 5T, 5N, 6T, 6N. The methods can be: WGBS, TWGBS, PBAT or SWIFT. The assembly has to be hg38.
#' @param wdir The output directory. If not provided, the working directory will be used. If non-existing, it will be created.
#' @export

benchmarking_plots <- function(benchmarking_data=NULL, user_data=NULL, verbose=TRUE, wdir=NULL, use.col=T){

  if (is.null(wdir)){
    if (verbose)
      cat("Using the default working directory as an output directory.  \n")
    wdir <- getwd()
  } else if (dir.exists(wdir)){
    if (verbose)
    cat("Using the", wdir, "as an output directory.  \n")
  } else {
    dir.create(file.path(wdir))
    cat("Directory", wdir, "doesn't exist. Creating it and using it as an output directory.   \n")
  }

  if (is.null(benchmarking_data)){
    if (verbose)
      cat("No benchmarking data vas provided. Using the package's own. \n")
    benchmarking_data <- bsmadez::benchmarking_data
  }

  df <- bsmadez::Preprocess.MethData(benchmarking_data=benchmarking_data, user_data = user_data)

  methods <- names(df)
  for (method in methods){
    df[[method]]$pipeline <- factor(df[[method]]$pipeline)
    levels <-  c(
      "methCT-OTP" = "OTP",
      "methCT-DKFZ_B060" = "DKFZ.B060",
      "methCT-BioMedX" = "BioMedX",
      "gemBS-IHEC" = "IHEC",
      "bistro" = "bistro",
      "GSNAP-UdS.NT" = "UdS.NT",
      "GSNAP-UdS.TT" = "UdS.TT",
      "LaJolla-WGBS" = "UCSD",
      "DKFZ_B370" = "DKFZ.C010",
      "BAT" = "FLI",
      "gemBS-2" = "cnag_gembs",
      "Bismark" = "bismark")

    df[[method]]$pipeline <- forcats::fct_recode(df[[method]]$pipeline, !!!levels)

  bsmadez::epiGB.MethPlot(epiGB.data = df[[method]], pdf = T,
                          pdf_path = wdir, filename = paste0(method, "_MethPlot"), use.col=use.col)

  bsmadez::epiGB.DistPlot(epiGB.data = df[[method]], pdf = T,
                          pdf_path = wdir, filename = paste0(method, "_DistPlot"), use.col=use.col, width = 7, height = 8)
  #bsmadez::epiGB.CorrPlot(epiGB.data = df[[method]], pdf = T, pdf_path = wdir, filename = paste0(method, "_CorrPlot"))
  bsmadez::epiGB.Heatmap(epiGB.data = df[[method]], pdf = T,
                         pdf_path = wdir, filename =paste0(method, "_Heatmap"))
  }

}




#
# ######F####unctions############
# wdir <- "Z:/Reka/30_EpiGB/results/plots_finalized/"
# files <- dir(path="Z:/Reka/30_EpiGB/results/RScripts.Lilian/", full.names = T)
# files <- files[-grep("Master_File", files)]
#
# for (file in files){
#   source(file)
# }
#
# methods <- c("WGBS", "TWGBS")
#
# # Files to read in:
#
# epiGB.df <- readRDS(file = "V:/Reka/30_EpiGB/Data/gold_standard/benchmarking_data_new.rds")
#
# region_annotation <- readRDS(file = "V:/Reka/30_EpiGB/Data/gold_standard/region_annotation.rds")
#
# consensus_calls <- readRDS("V:/Reka/30_EpiGB/Data/gold_standard/gold_standard_calls_means.rds")
#
# # add options for more methods
# df <- Preprocess.MethData(epiGB.df = epiGB.df, consensus.df = consensus_calls)
#
# for (method in methods){
# epiGB.MethPlot(epiGB.data = df[[method]], pdf = T, pdf_path = wdir, filename = paste0(method, "_MethPlot"))
# epiGB.DistPlot(epiGB.data = df[[method]], pdf = T, pdf_path = wdir, filename = paste0(method, "_DistPlot"))
# epiGB.CorrPlot2(epiGB.data = df[[method]], pdf = T, pdf_path = wdir, filename = paste0(method, "_CorrPlot"))
# epiGB.Heatmap(epiGB.data = df[[method]], pdf = T, pdf_path = wdir, filename =paste0(method, "_Heatmap"))
# }
#

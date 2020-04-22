
#' @title Preprocessing Methrix Data
#'
#' @description The following function is used to preprocess whole-genome bisulfite sequencing and
#' tagmented wgbs
#' benchmarking and optionally submitted user data for subsequent plotting.
#' @param epiGB.df A methrix object containing the bechmarking data.
#' @param user_data Optional user data input in a methrix object. The sample annotation should have the following columns: sample, method, assembly and pipeline.
#' The samples available are: 5T, 5N, 6T, 6N. The methods can be: WGBS, TWGBS, PBAT or SWIFT. The assembly has to be hg38.
#' @return A list with datasets processed for plotting. One list element contains the methylation calls and coverage information for one method.
#' @import dplyr
#' @import methrix
#' @export


Preprocess.MethData <- function(benchmarking_data = NULL,  user_data = NULL, verbose=TRUE){

  if (is.null(benchmarking_data)){
    benchmarking_data <- bsmadez::benchmarking_data
    if (verbose)
    cat("The benchmarking data was loaded. \n")
  }
  if(class(benchmarking_data)!="methrix"){
    stop("The benchmarking data has to be a methrix object. \n")
  }

  if (is.null(user_data)){
    cat("No user data was provided. \n")
  }
  if (verbose)
    cat("Checking the validity of the data.... \n")

  if (!(is.null(user_data))){
    if(class(user_data)!="methrix"){
      stop("The user data has to be a methrix object. \n")
    }

    if (!all(colnames(user_data@colData) %in% c("sample","method","assembly","pipeline"))){
      stop("The sample annotation of the user data has to have the following columns: sample, method, assembly, pipeline.  \n")
    }
    if (!all(as.character(unique(user_data@colData$sample)) %in% c("5N", "5T", "6N", "6T"))){
      stop("Unknown samples in the user data.  \n")
    }
    if (!all(unique(user_data@colData$method) %in% c("WGBS", "TWGBS", "PBAT", "SWIFT"))){
      stop("Unknown methods in the user  data. \n ")
    }

    if (!all(unique(user_data@colData$assembly) %in% c("hg38") )){
      stop("The assembly of the user data has to be hg38. \n")
    }

  }

  #sanity checks for the BM data
  if (!all(as.character(unique(benchmarking_data@colData$sample)) %in% c("5N", "5T", "6N", "6T"))){
    stop("Unknown samples in the bechmarking data.  \n")
  }

  if (!all(c("5N", "5T", "6N", "6T") %in% unique(benchmarking_data@colData$sample))){
    stop("Not all the samples are present in the benchmarking dataset.  \n")
  }

  if (!all(unique(benchmarking_data@colData$method) %in% c("WGBS", "TWGBS", "PBAT", "SWIFT"))){
    stop("Unknown methods in the bechmarking  data.  \n")
  }

  if (!all(unique(benchmarking_data@colData$assembly) %in% c("hg38"))){
    stop("The assembly of the bechmarking data has to be hg38. \n")
  }
  if (!(is.null(user_data))){
    if (any(user_data@colData$pipeline %in% benchmarking_data@colData$pipeline)){

      stop("One or more of the pipelines are already in the original benchmarking data. Please rename them. \n")
    }
}
if (verbose)
  cat("Checking is done, everything seems good. \n")

  region_annotation <- bsmadez::region_annotation
regions <- region_annotation[,c("seqnames", "probe_pos")]
regions$end <- regions$probe_pos+1
colnames(regions) <-c("chr", "start", "end")
regions <- GenomicRanges::makeGRangesFromDataFrame(regions, keep.extra.columns = T, ignore.strand = T)

if (!is.null(user_data)){
user_data <- methrix::subset_methrix(user_data, regions = regions)

if (nrow(user_data@elementMetadata)!=47){
  stop("Not all the regions were present in the user data. \n")
}

if (verbose)
  cat("Combining benhcmarking data and user data. \n")
###change it to methrix cbind
benchmarking_data <- SummarizedExperiment::cbind(benchmarking_data, user_data)

}


beta_mat <- assays(benchmarking_data)[["beta"]]
reg <- benchmarking_data@elementMetadata
reg <- as.data.frame(reg[,c("chr", "start")])

reg <- dplyr::inner_join(reg, region_annotation[,c("seqnames", "probe_pos", "locus_identifier")], by=c("start"="probe_pos", "chr"="seqnames"))

colnames(beta_mat)<- rownames(as.data.frame(benchmarking_data@colData))
colnames(beta_mat) <- paste0(colnames(beta_mat), "_meth")
reg <- cbind(reg, beta_mat)
cov_mat <- assays(benchmarking_data)[["cov"]]
colnames(cov_mat)<- rownames(as.data.frame(benchmarking_data@colData))
colnames(cov_mat) <- paste0(colnames(cov_mat), "_cov")
reg <- cbind(reg, cov_mat)
reg <- epiGB.melt(reg)

consensus.df <- bsmadez::consensus_calls
consensus.df <- data.table::rbindlist(consensus.df, T, idcol = "samples")
consensus.df$samples <- gsub("(X)([[:alnum:]]{2})(_meth)", "\\2", consensus.df$samples)



epiGB.df <- dplyr::left_join(x = reg, y = consensus.df[,c("samples", "locus_identifier", "lower", "upper", "mean")],
                            by=c("sample"="samples", "locus_identifier"="locus_identifier"))

# Splitting the dataset into two based on method:
# (Add option for user_data if it contains more than one method)
# Add user data to the list?:
# epiGB.df <- list(WGBS = epiGB.df[epiGB.df$method == "WGBS",], TWGBS = epiGB.df[epiGB.df$method == "TWGBS",], user_data)
# order locus_identifier numerically:

epiGB.df$locus_identifier <- factor(epiGB.df$locus_identifier,
                                    levels = c(paste0("mandatory_", seq(1:16)),
                                               paste0("recommended_", seq(1:32))))

methods <- unique(epiGB.df$method)

epiGB.df <- lapply(methods , function(x) epiGB.df[epiGB.df$method==x,])
names(epiGB.df) <- methods

epiGB.df <- lapply(X = epiGB.df, FUN = epiGB.dist)


return(epiGB.df)

}

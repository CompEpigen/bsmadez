
#' @title Data Melting
#'
#' @description The following function is used to melt the preprocessed benchmark and
#' optionally provided user data to a format that can be used by ggplot.
#' @param value The preprocessed data from submitted by epiGB.DataStructure.
#' @export
#' @import reshape2
#' @import dplyr


epiGB.melt <- function(df=NULL){
  tb_excluded <- colnames(df)[!colnames(df) %in%
                                   c("chr", "start", "locus_identifier")]
  meth_excluded <- tb_excluded[-grep("_meth", tb_excluded)]
  cov_excluded <- tb_excluded[-grep("_cov", tb_excluded)]

  #loop over the columns we want to change and make them numerical
  for (i in tb_excluded){
    df[,i]<- as.numeric(as.character(df[,i]))
  }

  #two step melting to melt both coverage and methylation
  mdf_1 <- reshape2::melt(df, id.vars = c("chr", "start", "locus_identifier"), measure.vars = cov_excluded)
  mdf_2 <- reshape2::melt(df, id.vars = c("chr", "start", "locus_identifier"), measure.vars = meth_excluded)


  #the "_meth" and "_cov" is subsituted before merge so the merge can happen successfully
  mdf_1$variable <- gsub("_meth", "", mdf_1$variable)
  mdf_2$variable <- gsub("_cov", "", mdf_2$variable)


  names_regexpr <- "hg38.([T]?WGBS).([5-6][N-T]).([[:alnum:]]+.?[[:alnum:]]*)"
  mdf_1$pipeline <- gsub(names_regexpr, "\\3", mdf_1$variable)
  mdf_1$method <- gsub(names_regexpr, "\\1", mdf_1$variable)
  mdf_1$sample <- gsub(names_regexpr, "\\2", mdf_1$variable)
  mdf_2$pipeline <- gsub(names_regexpr, "\\3", mdf_2$variable)
  mdf_2$method <- gsub(names_regexpr, "\\1", mdf_2$variable)
  mdf_2$sample <- gsub(names_regexpr, "\\2", mdf_2$variable)

  epiGB.df <- dplyr::inner_join(mdf_1, mdf_2, by = c("chr", "start", "locus_identifier", "variable", "pipeline", "method", "sample"))


  epiGB.df <- epiGB.df[c("chr", "start", "locus_identifier",  "method", "pipeline",  "sample", "value.x", "value.y")]
  colnames(epiGB.df) <- c("chr", "start", "locus_identifier", "method", "pipeline", "sample", "Methylation", "Coverage")
  #melted_df <- melted_df[, c("sample")]

  #epiGB.df <- merge(melted_df, mgv, by = c("sample", "regions"))
  #epiGB.df <<- epiGB.df
  return(epiGB.df)
}

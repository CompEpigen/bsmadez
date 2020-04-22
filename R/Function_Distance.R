
#' @title Distance Calculation
#'
#' @description The following function is used to calculate the absolute and
#' relative distance of benchmark and optionally provided user data from the consensus corridor.
#' @param df The preprocessed and melted benchmark and optionally provided user data.
#' @return Returns with a data frame with calculated absolute and relative distances
#' @export


epiGB.dist <- function(df=NULL){


  df$range <- with(df, df$lower <= df$Methylation &
                     df$upper >= df$Methylation)

  for (i in (1:nrow(df))){
    if (is.na(df[i, "range"])){
      df[i, "abs_dist"] <- NA
      df[i, "rel_dist"] <- NA
    } else if (df[i, "range"] == TRUE){
      df[i, "abs_dist"] <- 0
      df[i, "rel_dist"] <- 0
    }else{
      m = abs(df$upper[i] - df$Methylation[i])
      n = abs(df$lower[i] - df$Methylation[i])
      if ((m < n) == TRUE){
        df[i, "abs_dist"] <- m
        df[i, "rel_dist"] <- (df$upper[i] - df$Methylation[i])
      }else{
        df[i, "abs_dist"] <- n
        df[i, "rel_dist"] <- (df$lower[i] - df$Methylation[i])
      }
    }
  }
  return(df)
}

discrepancy_score <- function(m, comp_type=c("absolute", "relative"), value=c("mean", "median"),
                              file=NULL, n_threads=4, sample="all", preprocessed=F, folder =NULL, precalc_only=F, pipelines_o=NULL){


  pipelines <- unique(paste(m@colData$method, m@colData$pipeline, sep="&"))



  if (is.null(pipelines_o)){
  pipelines <- pipelines[order(pipelines)]} else if (all(pipelines_o %in% pipelines)){
    pipelines <- pipelines_o
  }

  #browser()
  combinations <- combn(pipelines, 2)
  if (comp_type=="absolute"){
    final_data <- as.data.frame(matrix(NA, nrow = length(pipelines)^2, ncol = 3))
    colnames(final_data) <- c("pipeline1", "pipeline2", "diff")
  } else if (comp_type=="relative"){
    final_data <- as.data.frame(matrix(NA, nrow = length(pipelines)^2, ncol = 4))
    colnames(final_data) <- c("pipeline1", "pipeline2", "diff_pos", "diff_neg")
  } else {
    stop("Wrong comparison type")
  }
  final_data$pipeline1 <- rep(pipelines, each=length(pipelines))
  final_data$pipeline2 <- rep(pipelines, length(pipelines))

  #library(doParallel)
  #cl <- makeCluster(n_threads)
  #registerDoParallel(cl)

  #foreach(ccs=as.list(as.data.frame(combinations))) %dopar% {
  for (ccs in as.list(as.data.frame(combinations))){

    if (preprocessed & !is.null(folder)){

      if (comp_type=="relative"){
        read_in_1 <- paste0(folder, "_", as.character(ccs[1]), "_", as.character(ccs[2]), "_relative_", "pos", "_", value,".txt")
        read_in_2 <- paste0(folder, "_", as.character(ccs[1]), "_", as.character(ccs[2]), "_relative_", "neg", "_", value,".txt")
        ####other order, one has to change the directions!
        read_in_12 <- paste0(folder, "_", as.character(ccs[2]), "_", as.character(ccs[1]), "_relative_", "pos", "_", value,".txt")
        read_in_22 <- paste0(folder, "_", as.character(ccs[2]), "_", as.character(ccs[1]), "_relative_", "neg", "_", value,".txt")

        if (file.exists(read_in_1) & file.exists(read_in_2)){
          pp1 <- read.table(file = read_in_1, header = T)
          pp2 <- read.table(file = read_in_2, header = T)

          if (sample=="all"){
            pp1 <- sum(pp1)/4
            pp2 <- sum(pp2)/4
            pp <- list(pp1, pp2)
            names(pp) <- c("pos", "neg")

          } else if (sample %in% samples){
            pp1 <- pp1[1, paste0("X", sample)]
            pp2 <- pp2[1, paste0("X", sample)]
            pp <- list(pp1, pp2)
            names(pp) <- c("pos", "neg")
          } else {
            stop("Wrong sample name")
          }
        }
        else if (file.exists(read_in_12) & file.exists(read_in_22)){

          pp1 <- read.table(file = read_in_12, header = T)
          pp2 <- read.table(file = read_in_22, header = T)
        #  pp1 <- -pp1
        #  pp2 <- -pp2

           if (sample=="all"){
             pp1 <- sum(pp1)/4
             pp2 <- sum(pp2)/4
             pp <- list(pp1, pp2)
             names(pp) <- c("pos", "neg")

           } else if (sample %in% samples){
             pp1 <- pp1[1, paste0("X", sample)]
             pp2 <- pp2[1, paste0("X", sample)]
             pp <- list(pp1, pp2)
             names(pp) <- c("pos", "neg")
           } else {
             stop("Wrong sample name")

           }
      } else {
          pp <- pairwise_discrepancy(m, pipeline1 = as.character(ccs[1]), pipeline2 = as.character(ccs[2]),
                                     comp_type=comp_type, sample=sample, value=value, file=file, precalc_only = precalc_only)
        }
      } else if (comp_type=="absolute"){
        read_in <- paste0(folder, "_", as.character(ccs[1]), "_", as.character(ccs[2]), "_absolute_", value,".txt")
        read_in2 <- paste0(folder, "_", as.character(ccs[2]), "_", as.character(ccs[1]), "_absolute_", value,".txt")
        #read_in <- paste0(folder, "_", as.character(ccs[1]), "_", as.character(ccs[2]), "_absolute_",".txt")
        if (file.exists(read_in) | file.exists(read_in2)){
         # | file.exists(read_in2)){
          if (file.exists(read_in))
            pp <- read.table(file = read_in, header = T)
            else
              pp <- read.table(file = read_in2, header = T)
          #pp <- ifelse(file.exists(read_in), read.table(file = read_in, header = T), read.table(file = read_in2, header = T))

          if (sample=="all"){
            pp <- sum(pp)/4
          } else if (sample %in% samples){
            pp <- pp[1, paste0("X", sample)]
          } else {
            stop("Wrong sample name")
          }
        } else {
          pp <- pairwise_discrepancy(m, pipeline1 = as.character(ccs[1]), pipeline2 = as.character(ccs[2]),
                                     comp_type=comp_type, sample=sample, value=value, file=file, precalc_only = precalc_only)}
      }
    } else {
      pp <- pairwise_discrepancy(m, pipeline1 = as.character(ccs[1]), pipeline2 = as.character(ccs[2]),
                                 comp_type=comp_type, sample=sample, value=value, file=file, precalc_only = precalc_only)
    }
    #final_data$sample[final_data$pipeline1==ccs[1] & final_data$pipeline2==ccs[2]] <- names(pp)
    if (comp_type=="absolute"){

      final_data$diff[final_data$pipeline1==as.character(ccs[1]) & final_data$pipeline2==as.character(ccs[2])] <- pp

    } else if (comp_type=="relative"){

      final_data$diff_pos[final_data$pipeline1==as.character(ccs[1]) & final_data$pipeline2==as.character(ccs[2])] <- pp$pos
      final_data$diff_neg[final_data$pipeline2==as.character(ccs[1]) & final_data$pipeline1==as.character(ccs[2])] <- pp$neg
    }
    print(ccs)
  }
  if (comp_type=="absolute"){
    #final_data$diff <- 1-final_data$diff
    final_data$diff[final_data$pipeline1==final_data$pipeline2] <- 0
  } else if (comp_type=="relative"){
    #final_data$diff_pos <- 1-final_data$diff_pos
    #final_data$diff_neg <- (-1)-final_data$diff_neg


    #organize data for plotting.
    final_data$diff_pos[final_data$pipeline1==final_data$pipeline2] <- 0
    final_data$diff_pos[is.na(final_data$diff_pos)] <- final_data$diff_neg[is.na(final_data$diff_pos)]
  }

  final_data$pipeline1 <- factor(final_data$pipeline1, levels = unique(final_data$pipeline1))
  final_data$pipeline2 <- factor(final_data$pipeline2, levels = unique(final_data$pipeline1))


  return(final_data)
}

pairwise_discrepancy <- function(m, pipeline1=character(), pipeline2= character(),
                                 comp_type=c("absolute", "relative"), sample="all", value=c("mean", "median"),
                                 file=NULL, precalc_only=F){

  pair1 <- strsplit(pipeline1, split="&")
  pair2 <- strsplit(pipeline2, split="&")

  pair1_mat <- methrix::get_matrix(m = m[,which(m@colData$method==pair1[[1]][1] & m@colData$pipeline==pair1[[1]][2])])
  pair2_mat <- methrix::get_matrix(m = m[,which(m@colData$method==pair2[[1]][1] & m@colData$pipeline==pair2[[1]][2])])

  samples <- unique(m@colData$sample)
  #match(samples, gsub("(hg38).(T?WGBS).([56]{1}[NT]{1}).([[:alnum:]]+).?([[:alnum:]]*)",  "\\3", colnames(pair1_mat)))
  pair1_mat <- pair1_mat[,match(samples, gsub("(hg38).(T?WGBS).([56]{1}[NT]{1}).([[:alnum:]]+).?([[:alnum:]]*)",  "\\3", colnames(pair1_mat)))]
  pair2_mat <- pair2_mat[,match(samples, gsub("(hg38).(T?WGBS).([56]{1}[NT]{1}).([[:alnum:]]+).?([[:alnum:]]*)",  "\\3", colnames(pair2_mat)))]

  difference <- pair1_mat-pair2_mat
  colnames(difference) <-  paste0( pair1[[1]][2],"_",
                                   pair2[[1]][2],"_",
                                   gsub( "(hg38).(T?WGBS).([56]{1}[NT]{1}).([[:alnum:]]+).?([[:alnum:]]*)","\\3",colnames(pair1_mat)))

  #browser()
  if (!precalc_only & !is.null(file)){
    fwrite(as.data.table(difference), file = paste0(file, "_", pipeline1, "_", pipeline2, ".txt"), row.names = F, sep="\t")
  }

  if (comp_type=="absolute"){
    if (class(difference)=="DelayedMatrix"){
      if (value=="mean"){
        sum_diff <- DelayedMatrixStats::colMeans2(abs(difference), na.rm = T)
      } else {
        sum_diff <- DelayedMatrixStats::colMedians(abs(difference), na.rm = T)}
    } else {
      if (value=="mean"){
        sum_diff <- matrixStats::colMeans2(abs(difference), na.rm = T)
      } else {
        sum_diff <- matrixStats::colMedians(abs(difference), na.rm = T)}
    }
    #####maybe implement median instead of mean

    #sum_diff <- sum_diff/nrow(m)
    names(sum_diff) <- gsub("(hg38).(T?WGBS).([56]{1}[NT]{1}).([[:alnum:]]+).?([[:alnum:]]*)",  "\\3", colnames(pair1_mat))
    if (sample=="all"){
      if (!is.null(file)){
        write.table(t(as.data.frame(sum_diff)),
                    file = paste0(file, "_", pipeline1, "_", pipeline2, "_absolute",  "_", value,".txt"), row.names = F, col.names = T, sep = "\t")
      }
      sum_diff <- sum(sum_diff)/4

    } else if (sample %in% samples){
      sum_diff <- sum_diff[sample]
    } else {
      stop("Wrong sample name")
    }
    #
  }

  if (comp_type=="relative"){
    difference_pos <- difference
    difference_pos[difference_pos<0] <- NA
    difference_neg <- difference
    difference_neg[difference_neg>0] <- NA

    if (class(difference)=="DelayedMatrix"){
      if (value == "mean") {
        sum_diff_pos <- DelayedMatrixStats::colMeans2(difference_pos, na.rm = T)
        sum_diff_neg <- DelayedMatrixStats::colMeans2(difference_neg, na.rm = T)
      } else {
        sum_diff_pos <- DelayedMatrixStats::colMedians(difference_pos, na.rm = T)
        sum_diff_neg <- DelayedMatrixStats::colMedians(difference_neg, na.rm = T)
      }
    } else {
      if (value == "mean") {
        sum_diff_pos <- matrixStats::colMeans2(difference_pos, na.rm = T)
        sum_diff_neg <- matrixStats::colMeans2(difference_neg, na.rm = T)
      } else {
        sum_diff_pos <- matrixStats::colMedians(difference_pos, na.rm = T)
        sum_diff_neg <- matrixStats::colMedians(difference_neg, na.rm = T)
      }
    }
    names(sum_diff_pos) <- gsub("(hg38).(T?WGBS).([56]{1}[NT]{1}).([[:alnum:]]+).?([[:alnum:]]*)",  "\\3", colnames(pair1_mat))
    names(sum_diff_neg) <- gsub("(hg38).(T?WGBS).([56]{1}[NT]{1}).([[:alnum:]]+).?([[:alnum:]]*)",  "\\3", colnames(pair1_mat))
    #sum_diff <- list((sum_diff_pos/nrow(m)), (sum_diff_neg/nrow(m)))
    sum_diff <- list(sum_diff_pos, sum_diff_neg)
    names(sum_diff) <- c("pos", "neg")

    if (sample=="all"){
      if (!is.null(file)){
        for (direction in c("pos", "neg")){
          write.table(as.data.frame(t(sum_diff[[direction]])),
                      file = paste0(file, "_", pipeline1, "_", pipeline2, "_relative_", direction, "_", value,".txt"), row.names = F, col.names = T, sep = "\t")
        }
      }
      sum_diff <- lapply(sum_diff, function(x) sum(x)/4)
    } else if (sample %in% samples){
      sum_diff <- lapply(sum_diff, function(x) x[sample])
    } else {
      stop("Wrong sample name")
    }
  }

  return(sum_diff)
  gc()
}

clustering_disc <- function(m, discrepancy_score, comp_type=c("absolute", "relative"), what=c("pos", "neg"),
                            folder=NULL, value="mean", sample="all", file=NULL){

  if (comp_type=="absolute"){
    dist <- reshape2::dcast(discrepancy_score, pipeline1 ~ pipeline2 )
    rownames(dist) <- dist[,1]
    dist <- dist[,-1]

    dist <- dplyr::arrange(dist, -dplyr::row_number())[,ncol(dist):1]
    clusters <- hclust(as.dist(dist))


    score <- discrepancy_score(
        m,
        comp_type = comp_type,
        value = value,
        sample = sample,
        folder = folder,
        preprocessed = T,
        pipelines_o = clusters$labels[clusters$order]
      )
    return(score)
  } else if (comp_type=="relative"){
    #browser()
    if (what == "pos"){
      score <- discrepancy_score[,c("pipeline1", "pipeline2", "diff_pos")]
      score$diff_pos[score$diff_pos<0] <- NA
      colnames(score)[3] <- "diff"
    } else if (what == "neg"){
      score <- discrepancy_score[,c("pipeline1", "pipeline2", "diff_neg")]
      colnames(score)[3] <- "diff"
      score$diff <- -score$diff
    }

    dist <- reshape2::dcast(score, pipeline1 ~ pipeline2)
    rownames(dist) <- dist[,1]
    dist <- dist[,-1]
if (what == "pos"){
    dist <- dplyr::arrange(dist, -dplyr::row_number())[,ncol(dist):1]}
    clusters <- hclust(as.dist(dist))


    score_clus <- discrepancy_score(
      m,
      comp_type = comp_type,
      value = value,
      sample = sample,
      folder = folder,
      preprocessed = T,
      pipelines_o = clusters$labels[clusters$order],
      file=file,
     precalc_only = T
    )
    return(score_clus)
  }
}

discrepancy_plot <- function(discrepancy_score, comp_type=c("absolute", "relative"), name=""){
#browser()
discrepancy_score$pipeline1 <- gsub("&", "_", as.character(discrepancy_score$pipeline1))
discrepancy_score$pipeline2 <- gsub("&", "_", as.character(discrepancy_score$pipeline2))


  discrepancy_score$pipeline1 <- factor(discrepancy_score$pipeline1, levels = rev(unique(discrepancy_score$pipeline1)))
  discrepancy_score$pipeline2 <- factor(discrepancy_score$pipeline2, levels = rev(unique(discrepancy_score$pipeline2)))


  if (comp_type=="absolute"){


    g <- ggplot(discrepancy_score)+
      geom_tile(aes(pipeline1, pipeline2, fill=diff))+
      theme_bw()+
      theme(panel.border = element_blank(), legend.key = element_rect(fill = guide_colourbar(nbin = 10)), axis.text.x = element_text(angle=90, vjust=0.5))+
      scale_fill_gradientn(colors=c("#24325FFF", "#FAFD7CFF"), limits = c(0, 0.18), na.value = "white", breaks=seq(from = 0, to =0.3, by=0.05))+
      xlab("Pipelines")+ylab("Pipelines")+labs(title = paste0("Discrepancy plots in samples: ", name), fill="Difference")
  }

  if (comp_type=="relative"){

    g <- ggplot(discrepancy_score)+
      geom_tile(aes(pipeline1, pipeline2, fill=diff_pos))+
      theme_bw()+
      theme(panel.border = element_blank(), , axis.text.x = element_text(angle=90, vjust=0.5))+
      scale_fill_gradientn(colors=c("darkred", "white",  "darkblue"), limits=c(-0.2,0.2), breaks=c(seq(from = -0.4, to=0, by=0.05), 0, seq(from = 0, to = 0.4, by=0.05)))+
      xlab("Pipelines")+ylab("Pipelines")+labs(title = paste0("Discrepancy plots in samples: ", name), fill="Difference")

  }
  return(g)
}




discrepancy_manhattan <- function(m, pipeline1, pipeline2, comp_type=c("absolute", "relative"), genome="hg38"){

  if (genome=="hg38"){
    chr_lengths= c("chr1"=248956422,
                   "chr2"=242193529,
                   "chr3"=198295559,
                   "chr4"=190214555,
                   "chr5"=181538259,
                   "chr6"=170805979,
                   "chr7"=159345973,
                   "chr8"=145138636,
                   "chr9"=138394717,
                   "chr10"=133797422,
                   "chr11"=135086622,
                   "chr12"=133275309,
                   "chr13"=114364328,
                   "chr14"=107043718,
                   "chr15"=101991189,
                   "chr16"=90338345,
                   "chr17"=83257441,
                   "chr18"=80373285,
                   "chr19"=58617616,
                   "chr20"=64444167,
                   "chr21"=46709983,
                   "chr22"=50818468,
                   "chrX"=156040895,
                   "chrY"=57227415)
  } else {
    cat("The supplied genome is not supported yet. Please use hg38.")
  }

  adding_values <- cumsum(chr_lengths)
  adding_values <- data.table::shift(adding_values, n = 1L, fill=0, type="shift", give.names = T)
  names(adding_values) <- names(chr_lengths)
  adding_values["end"] <- cumsum(chr_lengths)[length(chr_lengths)]
  ####calculate difference for each position

  pair1 <- strsplit(pipeline1, split="_")
  pair2 <- strsplit(pipeline2, split="_")

  pair1_mat <- methrix::get_matrix(m = m[,which(m@colData$method==pair1[[1]][1] & m@colData$pipeline==pair1[[1]][2])])
  pair2_mat <- methrix::get_matrix(m = m[,which(m@colData$method==pair2[[1]][1] & m@colData$pipeline==pair2[[1]][2])])

  samples <- unique(m@colData$sample)
  #match(samples, gsub("(hg38).(T?WGBS).([56]{1}[NT]{1}).([[:alnum:]]+).?([[:alnum:]]*)",  "\\3", colnames(pair1_mat)))
  pair1_mat <- pair1_mat[,match(samples, gsub("(hg38).(T?WGBS).([56]{1}[NT]{1}).([[:alnum:]]+).?([[:alnum:]]*)",  "\\3", colnames(pair1_mat)))]
  pair2_mat <- pair2_mat[,match(samples, gsub("(hg38).(T?WGBS).([56]{1}[NT]{1}).([[:alnum:]]+).?([[:alnum:]]*)",  "\\3", colnames(pair2_mat)))]

  difference <- pair1_mat-pair2_mat

  positions <- m@elementMetadata

  positions$to_add <- 0
  for (chr in names(chr_lengths)){
    positions$to_add[positions$chr==chr] <- adding_values[chr]
  }
  positions <-setDT(as.data.frame(positions))

  positions <- positions[,position:=start+to_add]

  positions <- cbind(positions, as.data.frame(difference))
  center <- adding_values
  for (chr in names(adding_values)){
    center[chr] <- ((adding_values[which(names(adding_values)==chr)+1]-adding_values[chr])/2)+adding_values[chr]
  }

  val <- rep(c(0.8, 1), length(chr_lengths)/2)
  names(val) <- names(chr_lengths)
  browser()
  g <- list()


  g <- ggplot(positions) + geom_point(aes(position, hg38.WGBS.5N.DKFZ.C010, color=hg38.WGBS.5N.DKFZ.C010, alpha=chr), size=0.1)+theme_bw()+
    scale_alpha_manual(values = val)+
    theme(panel.border = element_blank(),
          panel.grid.major.x = ele,
          panel.grid.minor.x = element_blank(),
          legend.position = "none")+
    scale_x_continuous(labels=names(adding_values)[-length(adding_values)],
                       breaks=center[-length(adding_values)])+
    scale_color_gradientn(colours = c("darkred", "grey90", "darkblue"), breaks=c(-1, 0, 1))


  #png("Z:/Reka/30_EpiGB/results/plots/test_manhattan.png", units="px", width=4000, height=800, res=300)
  #g
  #dev.off()
  return(g)
}



colors <- ggsci::pal_rickandmorty(palette = "schwifty")(12)

levels <-  c(
  "methCT-OTP" = "OTP",
  "methCT-DKFZ_B060" = "DKFZ.B060",
  "methCT-DKFZ_B060" = "DKFZ-B060",
  "methCT-BioMedX" = "BioMedX",
  "gemBS-IHEC" = "IHEC",
  "Bistro" = "bistro",
  "GSNAP-UdS.NT" = "UdS.NT",
  "GSNAP-UdS.NT" = "UdS-NT",
  "GSNAP-UdS.TT" = "UdS.TT",
  "GSNAP-UdS.TT" = "UdS-TT",
  "LaJolla-WGBS" = "UCSD",
  "DKFZ_B370" = "DKFZ.C010",
  "DKFZ_B370" = "DKFZ-C010",
  "BAT" = "FLI",
  "gemBS-CNAG" = "cnag.gembs",
  "Bismark" = "bismark")

col_list <- c(
  "methCT-OTP" = colors[1],
  "methCT-DKFZ_B060" = colors[2],
  "methCT-BioMedX" = colors[3],
  "gemBS-IHEC" = colors[4],
  "Bistro" = colors[5],
  "GSNAP-UdS.NT" = colors[6],
  "GSNAP-UdS.TT" = colors[7],
  "LaJolla-WGBS" = colors[8],
  "DKFZ_B370" = colors[9],
  "BAT" = colors[10],
  "gemBS-2" = colors[11],
  "Bismark" = colors[12])

col_list2 <- c(
  "5"= colors[3],
  "6" = colors[5]
)

RES <- "V:/Reka/30_EpiGB/results/plots_20_01/"

m <- load_HDF5_methrix(dir ="N:/Internal/BS_benchmark/methylation_calls/read_in_hdf5/all_pipelines_16_01/")
m@colData$pipeline <- as.character(forcats::fct_recode(as.factor(m@colData$pipeline), !!!levels))

benchmarking_data <- bsmadez::benchmarking_data
benchmarking_data@colData$pipeline <- as.character(forcats::fct_recode(as.factor(benchmarking_data@colData$pipeline), !!!levels))

df <- bsmadez::Preprocess.MethData(benchmarking_data=benchmarking_data, user_data = NULL)
for (method in names(df)){
  df[[method]]$pipeline <- as.character(forcats::fct_recode(as.factor(df[[method]]$pipeline), !!!levels))
}

#
# ranking <- rank_twgbs %>% grank_twgbs <- rank.bySample(df)
# group_by(method, pipeline) %>% summarise(sum_dist=sum(sum_abs_dist))
#
#
#
# ranking$pipeline <- factor(ranking$pipeline,
#                            levels = unique(c(
#                              as.character(ranking$pipeline[ranking$method == "WGBS"][order(ranking$sum_dist[ranking$method ==
#                                                                                                               "WGBS"], decreasing = T)]),
#                              as.character(ranking$pipeline[ranking$method == "TWGBS"])
#                            )))
#
# ranking$pipeline <- factor(ranking$pipeline,
#                            levels = rev(levels(ranking$pipeline)))
#
# pdf("Y:/Internal/BS_benchmark/analysis/benchmarking_plots_14_07/Ranking.pdf", width = 7, height = 4)
# g <- ggplot()+
#   geom_point(data=ranking, aes(x = sum_dist, y = pipeline),
#                             size = 3*1.2, colour = "grey10", alpha = 0.7)+
#   geom_point(data=ranking, aes(x=sum_dist, y=pipeline, color=pipeline), size=3,  alpha = 0.7)+facet_grid(.~method)+
#   labs(title = "Ranking by method")+
#   labs(colour = "Pipeline")+
#   guides(colour = guide_legend(override.aes = list(alpha = 1)),
#          size = guide_legend(override.aes = list(alpha = 0.1)))+
#   xlab("Summarized distance")+
#   xlim(0, max(ranking$sum_dist)+1)+
#   ylab("Pipelines")+
#   theme_bw()+#this gets rid of the grey background
#   theme(strip.background =element_rect(fill="white"))+
#   theme(panel.grid.minor = element_blank(), legend.position = "none")
# if (!is.null(col_list)){
#   g <- g +  scale_color_manual(values=col_list, breaks=names(col_list))} else{
#     g <- g +  scale_color_brewer(palette="Spectral")
#   }
# print(g)
# dev.off()
#



mean_by_pipeline <- do.call(rbind.data.frame, df) %>%
  dplyr::group_by(method, pipeline, sample) %>%
  dplyr::summarize(av_dist = mean(abs_dist, na.rm=T), av_coverage = mean(Coverage, na.rm=T))

mean_by_pipeline$patient <- substr(mean_by_pipeline$sample, 1, 1)
mean_by_pipeline$type <- substr(mean_by_pipeline$sample, 2, 2)
mean_by_pipeline$type <- ifelse(mean_by_pipeline$type=="N", "Normal", "Tumor")
#pdf(paste0(RES, "/abs_distance.pdf"), width = 7, height = 4)

abs_dist <- summaryPlot_3(data=mean_by_pipeline, x = "pipeline", y = "av_dist", size = "av_coverage", group = "patient", types="type",
                          ylab="absolute distance", title="Avarage Absolute Distance from the Consensus Corridor", col_list=col_list2)
print(abs_dist)
#dev.off()
summaryPlot_3 <- function(data, x, y, size, group, types, ylab="", title="", col_list){

  g <- ggplot2::ggplot() + geom_point(data=data, aes(x = get(y), y = get(x),
                                                     size = get(size)*1.2, group = get(group), shape=get(types)),
                                      colour = "grey40", alpha = 0.7)+
    #,
    #position=position_dodge(width=0.5))+
    geom_point(data=data, aes(x = get(y), y = get(x),
                              colour = get(group), size = get(size), shape=get(types)),
               alpha = 0.7)+
    #, position=position_dodge(width=0.5))+
    facet_grid(.~method)+
    labs(title = title)+
    labs(colour = "Patient")+
    labs(size = "Coverage")+
    labs(shape="Sample Type")+
    guides(colour = guide_legend(override.aes = list(alpha = 1)),
           size = guide_legend(override.aes = list(alpha = 0.1)))+
    xlab(ylab)+
    xlim(0, max(data[,y]))+
    ylab("Pipelines")+
    theme_bw()+#this gets rid of the grey background
    #coord_flip()+#flip the axes
    theme(strip.background =element_rect(fill="white"))+
    theme(panel.grid.minor = element_blank())
  if (!is.null(col_list)){
    g <- g +  scale_color_manual(values=col_list, breaks=names(col_list))} else{
      g <- g +  scale_color_brewer(palette="Spectral")
    }

  return(g)
}


######Global methylation levels

selected_rows <- which(m@elementMetadata$chr %in% paste0("chr", 1:22))
m_2 <- m[selected_rows, ]
meth <- get_matrix(m_2, type="M")
means <- DelayedMatrixStats::colMeans2(meth, na.rm=T)

cov <- get_matrix(m_2, type="C")
coverage <- DelayedMatrixStats::colMeans2(cov, na.rm=T)

anno <- colData(m)
anno$means <- means
anno$coverage <- coverage
anno$patient <- substr( anno$sample, 1, 1)
anno$type <- substr( anno$sample, 2, 2)
anno$type <- ifelse(anno$type=="N", "Normal", "Tumor")
anno <- as.data.frame(anno)
anno$method <- factor(anno$method, levels=rev(levels(as.factor(anno$method))))

pdf(paste0(RES, "/global_methylation.pdf"), width = 7, height = 4)

g <- ggplot2::ggplot() + geom_point(data=anno, aes(x = means, y = pipeline,
shape=type, size=coverage*1.2),
  colour = "grey20", alpha = 0.7)+
#,
#position=position_dodge(width=0.5))+
geom_point(data=anno, aes(x = means, y = pipeline,
colour = patient,  shape=type, size=coverage),
alpha = 0.7)+
#, position=position_dodge(width=0.5))+
facet_grid(.~method)+
labs(title = "Mean methylation of autosomes")+
labs(colour = "Patient")+
labs(shape="Sample Type")+

guides(colour = guide_legend(override.aes = list(alpha = 1)),
size = guide_legend(override.aes = list(alpha = 0.1)))+
xlab("Mean")+
xlim(0.5, max(anno[,"means"]))+
ylab("Pipelines")+
theme_bw()+#this gets rid of the grey background
#coord_flip()+#flip the axes
theme(strip.background =element_rect(fill="white"))+
theme(panel.grid.minor = element_blank())+scale_color_manual(values=col_list2, breaks=names(col_list2))+
labs(size="Coverage")
print(g)

dev.off()

######mean ranking
methods <- c(names(df))
for (method in methods){

ranking <- df[[method]] %>% group_by(pipeline) %>% summarise(sum_dist=mean(abs_dist, na.rm=T))


ranking$pipeline <- factor(ranking$pipeline,
                           levels = unique(ranking$pipeline[order(ranking$sum_dist)]))
ranking$pipeline <- factor(ranking$pipeline,
                           levels = rev(levels(ranking$pipeline)))

pdf(paste0(RES, "/Ranking_mean_", method,  ".pdf"), width = 4, height = 4)
g <- ggplot()+
  geom_point(data=ranking, aes(x = sum_dist, y = pipeline),
             size = 3*1.2, colour = "grey10", alpha = 0.7)+
  geom_point(data=ranking, aes(x=sum_dist, y=pipeline, color=pipeline), size=3,  alpha = 0.7)+
  labs(title = method)+
  labs(colour = "Pipeline")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         size = guide_legend(override.aes = list(alpha = 0.1)))+
  xlab("Summarized distance")+
  xlim(0, max(ranking$sum_dist)+0.03)+
  ylab("Pipelines")+
  theme_bw()+#this gets rid of the grey background
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank(), legend.position = "none")
if (!is.null(col_list)){
  g <- g +  scale_color_manual(values=col_list, breaks=names(col_list))} else{
    g <- g +  scale_color_brewer(palette="Spectral")
  }
print(g)
dev.off()
}



for (method in methods){

  ranking <- df[[method]] %>% group_by(pipeline, sample) %>% summarise(sum_dist=mean(abs_dist, na.rm=T))
  sum_combined  <- ranking %>% group_by(pipeline) %>% summarise(sum_dist_comb=mean(sum_dist, na.rm=T))

  ranking$pipeline <- factor(ranking$pipeline,
                             levels = unique(sum_combined$pipeline[order(sum_combined$sum_dist_comb)]))
  ranking$pipeline <- factor(ranking$pipeline,
                             levels = rev(levels(ranking$pipeline)))

  ranking$patient <- substr(ranking$sample, 1, 1)
  ranking$type <- substr(ranking$sample, 2, 2)
  ranking$type <- ifelse(ranking$type=="N", "Normal", "Tumor")
  sum_combined  <- ranking %>% group_by(pipeline) %>% summarise(sum_dist=mean(sum_dist, na.rm=T))





  g <- ggplot2::ggplot() +
    stat_summary(data=ranking, aes(x = pipeline, y = sum_dist), fun.y="mean", geom="tile", width=0.6, height=0.001)+
    #geom_boxplot(data=ranking, aes(x=pipeline, y=sum_dist), color="grey70")+
  geom_point(data=ranking, aes(x = pipeline, y = sum_dist,
                                                     group = patient, shape=type),
                                      size = 3*1.2, colour = "grey40", alpha = 0.7)+
    #,
    #position=position_dodge(width=0.5))+
    geom_point(data=ranking, aes(x = pipeline, y = sum_dist, group=patient, shape=type, colour = patient), size=3,  alpha = 0.7)+
    #, position=position_dodge(width=0.5))+
    labs(title = method)+
    labs(colour = "Patient")+
    labs(size = "Coverage")+
    labs(shape="Sample Type")+
    guides(colour = guide_legend(override.aes = list(alpha = 1)),
           size = guide_legend(override.aes = list(alpha = 0.1)))+
    ylab("Avarage distance from the consensus corridor")+
    xlab("Pipelines")+
    ylim(0, max(ranking$sum_dist)+0.03)+
    theme_bw()+#this gets rid of the grey background
    coord_flip()+#flip the axes
    theme(strip.background =element_rect(fill="white"))+
    theme(panel.grid.minor = element_blank())
  if (!is.null(col_list)){
    g <- g +  scale_color_manual(values=col_list2, breaks=names(col_list2))} else{
      g <- g +  scale_color_brewer(palette="Spectral")
    }


  pdf(paste0(RES, "/Combined_ranking_", method,  ".pdf"), width = 6, height = 4)
  print(g)
  dev.off()

}



#####global coverage

n_covered_per_chr <- read.delim(paste0(RES, "/n_covered_per_chr.tsv"), stringsAsFactors=FALSE)
n_covered_per_chr <- n_covered_per_chr[n_covered_per_chr$chr %in% paste0("chr", 1:22),]
setDT(n_covered_per_chr)
global_summary <- n_covered_per_chr[,.(covered=sum(n_covered, na.rm=T), all=sum(total_CpGs, na.rm=T)), by=.(Sample_Name)]

name_regex <- "(hg[[:digit:]]{2})\\.(T?WGBS)\\.([56][NT])\\.(.*)"

global_summary$sample <- gsub(name_regex, "\\3", global_summary$Sample_Name)
global_summary$method <-  gsub(name_regex, "\\2", global_summary$Sample_Name)
global_summary$method <- factor(global_summary$method, levels=c("WGBS", "TWGBS"))
global_summary$assembly <-  gsub(name_regex, "\\1", global_summary$Sample_Name)
global_summary$pipeline <-  gsub(name_regex, "\\4", global_summary$Sample_Name)

global_summary$pipeline <- as.character(forcats::fct_recode(as.factor(global_summary$pipeline), !!!levels))
global_summary$patient <- substr(global_summary$sample, 1, 1)
global_summary$type <- substr(global_summary$sample, 2, 2)
global_summary$type <- ifelse(global_summary$type=="N", "Normal", "Tumor")
global_summary$covered <- global_summary$covered/1000000


pdf(paste0(RES, "/global_coverage.pdf"), width = 7, height = 4)

g <- ggplot2::ggplot() + geom_point(data=global_summary, aes(x = covered, y = pipeline,
                                                   shape=type), size=2*1.2,
                                    colour = "grey20", alpha = 0.7)+
  #,
  #position=position_dodge(width=0.5))+
  geom_point(data=global_summary, aes(x = covered, y = pipeline,
                            colour = patient,  shape=type), size=2,
             alpha = 0.7)+
  #, position=position_dodge(width=0.5))+
  facet_grid(.~method)+
  labs(title = "Global coverage")+
  labs(colour = "Patient")+
  labs(shape="Sample Type")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         size = guide_legend(override.aes = list(alpha = 0.1)))+
  xlab("CpGs covered (million)")+
  #xlim(0.5, max(anno[,"means"]))+
  ylab("Pipelines")+
  theme_bw()+#this gets rid of the grey background
  #coord_flip()+#flip the axes
  theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(size=8))+
  theme(panel.grid.minor = element_blank())+scale_color_manual(values=col_list2, breaks=names(col_list2))

print(g)

dev.off()


#####Depth of coverage
statistics_global <- read.delim(paste0(RES, "/global_MC_per_samp.tsv"))
statistics_global$sample <- gsub(name_regex, "\\3", statistics_global$Sample_Name)
statistics_global$method <-  gsub(name_regex, "\\2", statistics_global$Sample_Name)
statistics_global$method <- factor(statistics_global$method, levels=c("WGBS", "TWGBS"))
statistics_global$assembly <-  gsub(name_regex, "\\1", statistics_global$Sample_Name)
statistics_global$pipeline <-  gsub(name_regex, "\\4", statistics_global$Sample_Name)

statistics_global$pipeline <- as.character(forcats::fct_recode(as.factor(statistics_global$pipeline), !!!levels))
statistics_global$patient <- substr(statistics_global$sample, 1, 1)
statistics_global$type <- substr(statistics_global$sample, 2, 2)
statistics_global$type <- ifelse(statistics_global$type=="N", "Normal", "Tumor")
#global_summary$covered <- global_summary$covered/1000000

pdf(paste0(RES, "/coverage_depth.pdf"), width = 7, height = 4)


g <- ggplot2::ggplot() + geom_point(data=statistics_global, aes(x = mean_cov, y = pipeline,
                                                             shape=type), size=2*1.2,
                                    colour = "grey20", alpha = 0.7)+
  #,
  #position=position_dodge(width=0.5))+
  geom_point(data=statistics_global, aes(x = mean_cov, y = pipeline,
                                      colour = patient,  shape=type), size=2,
             alpha = 0.7)+
  #, position=position_dodge(width=0.5))+
  facet_grid(.~method, scales = "free")+
  labs(title = "Coverage depth")+
  labs(colour = "Patient")+
  labs(shape="Sample Type")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         size = guide_legend(override.aes = list(alpha = 0.1)))+
  xlab("Depth")+
  #xlim(0.5, max(anno[,"means"]))+
  ylab("Pipelines")+
  theme_bw()+#this gets rid of the grey background
  #coord_flip()+#flip the axes
  theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(size=8))+
  theme(panel.grid.minor = element_blank())+scale_color_manual(values=col_list2, breaks=names(col_list2))

print(g)

dev.off()

pdf(paste0(RES, "/global_meth.pdf"), width = 7, height = 4)


g <- ggplot2::ggplot() + geom_point(data=statistics_global, aes(x = mean_meth, y = pipeline,
                                                                shape=type), size=2*1.2,
                                    colour = "grey20", alpha = 0.7)+
  #,
  #position=position_dodge(width=0.5))+
  geom_point(data=statistics_global, aes(x = mean_meth, y = pipeline,
                                         colour = patient,  shape=type), size=2,
             alpha = 0.7)+
  #, position=position_dodge(width=0.5))+
  facet_grid(.~method)+
  xlim(0,1)+
  labs(title = "Methylation")+
  labs(colour = "Patient")+
  labs(shape="Sample Type")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         size = guide_legend(override.aes = list(alpha = 0.1)))+
  xlab("Global methylation")+
  #geom_vline(xintercept=global_meth[1], linetype="dashed")+
  #xlim(0.5, max(anno[,"means"]))+
  ylab("Pipelines")+
  theme_bw()+#this gets rid of the grey background
  #coord_flip()+#flip the axes
  theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(size=8))+
  theme(panel.grid.minor = element_blank())+scale_color_manual(values=col_list2, breaks=names(col_list2))

print(g)

dev.off()

global_meth <- c("5N"=77.61070238,
"5T"=75.52880953,
"6N"=77.41122619,
"6T"=76.6147381
)
global_meth <- global_meth/100

pdf(paste0(RES, "/global_meth_lines.pdf"), width = 7, height = 4)


g <- ggplot2::ggplot() + geom_point(data=statistics_global, aes(x = mean_meth, y = pipeline,
                                                                shape=type), size=2*1.2,
                                    colour = "grey20", alpha = 0.7)+
  #,
  #position=position_dodge(width=0.5))+
  geom_point(data=statistics_global, aes(x = mean_meth, y = pipeline,
                                         colour = patient,  shape=type), size=2,
             alpha = 0.7)+
  #, position=position_dodge(width=0.5))+
  facet_grid(.~method)+
  xlim(0.6,0.80)+
  labs(title = "Methylation")+
  labs(colour = "Patient")+
  labs(shape="Sample Type")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         size = guide_legend(override.aes = list(alpha = 0.1)))+
  xlab("Global methylation")+
  geom_vline(xintercept=global_meth, linetype="dashed")+
  #xlim(0.5, max(anno[,"means"]))+
  ylab("Pipelines")+
  theme_bw()+#this gets rid of the grey background
  #coord_flip()+#flip the axes
  theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(size=8))+
  theme(panel.grid.minor = element_blank())+scale_color_manual(values=col_list2, breaks=names(col_list2))

print(g)

dev.off()

statistics_global$difference <- NA
statistics_global$difference[statistics_global$sample=="5T"] <- statistics_global$mean_meth[statistics_global$sample=="5T"]-global_meth["5T"]
statistics_global$difference[statistics_global$sample=="6T"] <- statistics_global$mean_meth[statistics_global$sample=="6T"]-global_meth["6T"]

statistics_global$difference[statistics_global$sample=="5N"] <- statistics_global$mean_meth[statistics_global$sample=="5N"]- global_meth["5N"]
statistics_global$difference[statistics_global$sample=="6N"] <- statistics_global$mean_meth[statistics_global$sample=="6N"]-global_meth["6N"]


      st_gl_WGBS <- statistics_global[statistics_global$method=="WGBS",]
st_gl_TWGBS <- statistics_global[statistics_global$method=="TWGBS",]

st_gl_WGBS$pipeline <- factor(st_gl_WGBS$pipeline,
levels=st_gl_WGBS$pipeline[st_gl_WGBS$sample=="5T"][order(st_gl_WGBS$difference[st_gl_WGBS$sample=="5T"])])

st_gl_TWGBS$pipeline <- factor(st_gl_TWGBS$pipeline,
                              levels=st_gl_TWGBS$pipeline[st_gl_TWGBS$sample=="5T"][order(st_gl_TWGBS$difference[st_gl_TWGBS$sample=="5T"])])


pdf(paste0(RES, "/global_meth_diff_WGBS.pdf"), width = 7, height = 4)


g <- ggplot2::ggplot() + geom_point(data=st_gl_WGBS, aes(x = difference, y = pipeline,
                                                                shape=type), size=2*1.2,
                                    colour = "grey20", alpha = 0.7)+
  #,
  #position=position_dodge(width=0.5))+
  geom_point(data=st_gl_WGBS, aes(x = difference, y = pipeline,
                                         colour = patient,  shape=type), size=2,
             alpha = 0.7)+
  #, position=position_dodge(width=0.5))+
  #facet_grid(.~method)+
  #xlim(0.6,0.80)+
  labs(title = "WGBS")+
  labs(colour = "Patient")+
  labs(shape="Sample Type")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         size = guide_legend(override.aes = list(alpha = 0.1)))+
  xlab("Error of global methylation measurement")+
  geom_vline(xintercept=0, linetype="dashed")+
  #xlim(0.5, max(anno[,"means"]))+
  ylab("Pipelines")+
  theme_bw()+#this gets rid of the grey background
  #coord_flip()+#flip the axes
  theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(size=12))+
  theme(panel.grid.minor = element_blank())+scale_color_manual(values=col_list2, breaks=names(col_list2))

print(g)

dev.off()

pdf(paste0(RES, "/global_meth_diff_TWGBS.pdf"), width = 7, height = 4)


g <- ggplot2::ggplot() + geom_point(data=st_gl_TWGBS, aes(x = difference, y = pipeline,
                                                         shape=type), size=2*1.2,
                                    colour = "grey20", alpha = 0.7)+
  #,
  #position=position_dodge(width=0.5))+
  geom_point(data=st_gl_TWGBS, aes(x = difference, y = pipeline,
                                  colour = patient,  shape=type), size=2,
             alpha = 0.7)+
  #, position=position_dodge(width=0.5))+
  #facet_grid(.~method)+
  #xlim(0.6,0.80)+
  labs(title = "TWGBS")+
  labs(colour = "Patient")+
  labs(shape="Sample Type")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         size = guide_legend(override.aes = list(alpha = 0.1)))+
  xlab("Error of global methylation measurement ")+
  geom_vline(xintercept=0, linetype="dashed")+
  #xlim(0.5, max(anno[,"means"]))+
  ylab("Pipelines")+
  theme_bw()+#this gets rid of the grey background
  #coord_flip()+#flip the axes
  theme(strip.background =element_rect(fill="white"), axis.text.x = element_text(size=12))+
  theme(panel.grid.minor = element_blank())+scale_color_manual(values=col_list2, breaks=names(col_list2))

print(g)

dev.off()


score <- discrepancy_score(m, comp_type = "absolute", value="mean", folder = "N:/Internal/BS_benchmark/analysis/discrepancy_scores/discrepancy_score_abs_mean", preprocessed = T)
score_cluster <- clustering_disc(m, score, comp_type = "absolute", what = "pos", folder = "N:/Internal/BS_benchmark/analysis/discrepancy_scores/discrepancy_score_abs_mean")
g <- discrepancy_plot(discrepancy_score = score_cluster, comp_type = "absolute", name = "discrepancy_score_22_01")

ggsave(g, filename = paste0(RES, "discrepancy_plot_abs.pdf"))







######################
#small test
 m <- m[1:100000,]
#
m@assays[[1]][m@assays[[2]]==0] <- as.double(NA)
m@assays[[2]][m@assays[[2]]==0] <- as.integer(NA)

statistics <- get_stats(m, per_chr = F)

name_regex <- "(hg[[:digit:]]{2})\\.(T?WGBS)\\.([56][NT])\\.(.*)"
statistics$sample <- gsub(name_regex, "\\3", statistics$Sample_Name)
statistics$method <-  gsub(name_regex, "\\2", statistics$Sample_Name)
statistics$method <- factor(statistics$method, levels=c("WGBS", "TWGBS"))
statistics$assembly <-  gsub(name_regex, "\\1", statistics$Sample_Name)
statistics$pipeline <-  gsub(name_regex, "\\4", statistics$Sample_Name)

statistics$pipeline <- as.character(forcats::fct_recode(as.factor(statistics$pipeline), !!!levels))
statistics$patient <- substr(statistics$sample, 1, 1)
statistics$type <- substr(statistics$sample, 2, 2)
statistics$type <- ifelse(statistics$type=="N", "Normal", "Tumor")
##############

#DARs with CSAW

# Differentially accsesible regions from BAM files

library(csaw)
library(edgeR)
library(rtracklayer)
library(GenomicRanges)
library(org.At.tair.db)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(biomaRt)
library(rtracklayer)
library(ChIPseeker)
library(clusterProfiler)



############# VARIABLES TO ADJUST AS NEEDED ############################################################################

# Input files and names

experiment <- "noPwt_v_noPmut" # ID for output files

bam.files <- c(
  "wt_noP_1.nuclear.bam",
  "wt_noP_2.nuclear.bam",
  "mut_noP_1.nuclear.bam",
  "mut_noP_2.nuclear.bam"
)

condiciones <- c(rep("wt_noP", 2),
                 rep("mut_noP", 2)) # conditions for edgeR matrix (to work with duplicates)


############# ACTUAL CODE  ############################################################################################

# csaw general parameters

pe.param <-
  readParam(max.frag = 1000, pe = "both") # general parameters

# window counting
ventanasCounts <-
  windowCounts(bam.files = bam.files,
               width = 75,
               param = pe.param)


# Global enrichment filtering
bin.size <- 500
binned <-
  windowCounts(
    bam.files = bam.files,
    bin = T,
    width = bin.size,
    param = pe.param
  )


filter.stat <-
  filterWindowsGlobal(ventanasCounts, background = binned)
keep <- filter.stat$filter > log2(2)
sum(keep)

filtered.data <- ventanasCounts[keep,]


pdf(file = paste0(experiment,"_background_distribution.pdf"), width = 6, height = 6) 

par(mfrow = c(1, 1))
hist(
  filter.stat$back.abundances,
  xlab = "Adjusted bin log-CPM",
  breaks = 100,
  main = "",
  col = "grey80"
)
global.bg <- filter.stat$abundances - filter.stat$filter
abline(v = global.bg[1], col = "red", lwd = 2)
abline(v = global.bg[1] + log2(2),
       col = "blue",
       lwd = 2)

legend(
  "topright",
  lwd = 2,
  col = c('red', 'blue'),
  legend = c("Background", "log2FC = 2")
)

dev.off()

# Normalizing by TMM

filtered.data <- normFactors(filtered.data, se.out = filtered.data)

effORI  <- filtered.data$norm.factors
compORI <- normFactors(binned, se.out = F)
adjcORI <- cpm(asDGEList(binned), log = TRUE)

pdf(file = paste0(experiment,"_samplesExploration.pdf"), width = 6, height = 6) 

par(mfrow = c(1, 3), mar = c(5, 4, 2, 1.5))

for (i in seq_len(length(bam.files) - 1)) {
  adjc <- adjcORI[, c(1, 1 + i)]
  eff <- effORI[c(1, 1 + i)]
  comp <- compORI[c(1, 1 + i)]
  smoothScatter(
    x = rowMeans(adjc),
    y = adjc[, 1] - adjc[, 2],
    xlab = "A",
    ylab = "M",
    main = paste("1 vs", i + 1)
  )
  abline(h = log2(eff[1] / eff[2]), col = "red")
  abline(h = log2(comp[1] / comp[2]),
         col = "red",
         lty = 2)
}

dev.off()

# EDGER #####

y <- asDGEList(filtered.data)

# matrix design
design <- model.matrix(~ factor(condiciones))
colnames(design) <- c("intercept", "condiciones")
design

# dispersion estimation

y <- estimateDisp(y, design)
summary(y$trended.dispersion)

fit <- glmQLFit(y, design, robust = T)
summary(fit$var.post)

pdf(file = paste0(experiment,"_edgeRmodelExploration.pdf"), width = 6, height = 6)
par(mfrow = c(1, 2))
o <- order(y$AveLogCPM)
plot(
  y$AveLogCPM[o],
  sqrt(y$trended.dispersion[o]),
  type = "l",
  lwd = 2,
  ylim = c(0, 1),
  xlab = expression("Ave." ~ Log[2] ~ "CPM"),
  ylab = ("Biological coefficient of variation")
)
plotQLDisp(fit)
dev.off()

# Testing for DB windows #####

results <- glmQLFTest(fit, contrast = c(0, 1))

head(results$table)

rowData(filtered.data) <-
  cbind(rowData(filtered.data), results$table)

# MDS Plots

sampleidx <- c("_1", "_2")
Labels <- paste0(condiciones, sampleidx)


pdf(file = paste0(experiment,"_MDSplots.pdf"), width = 6, height = 6)
par(mfrow = c(2, 2), mar = c(5, 4, 2, 2))
adj.counts <- cpm(y, log = TRUE)
for (top in c(100, 500, 1000, 5000)) {
  out <-
    plotMDS(
      adj.counts,
      main = top,
      col = c("blue", "blue", "darkorange", "darkorange"),
      labels = Labels,
      top = top
    )
}
dev.off()

# Correction for multiple testing

merged_results <-
  mergeResults(
    filtered.data,
    results$table,
    tol = 150,
    merge.args = list(max.width = 5000)
  )


# log2FC and FDR filtering of the results and export ####

# FDR filtering
tf1 <- merged_results
tf1 <- tf1[tf1$combined$FDR < 0.05,]

# Removal of uncertain directionality differential windows
table(tf1$combined$direction)
tf1 <- tf1[tf1$combined$direction != "mixed",]

# log2FC filtering

temp1ups <- tf1[tf1$best$logFC > 0.8,]
temp1dws <- tf1[tf1$best$logFC < -0.8,]

ups <- temp1ups$regions
ups$score <- -10 * log10(temp1ups$combined$FDR)
names(ups) <- paste0("UpDar", 1:length(ups))
upsbedname <- paste0("upDARS_", experiment, ".bed")
export(ups, upsbedname)

downs <- temp1dws$regions
downs$score <- -10 * log10(temp1dws$combined$FDR)
names(downs) <- paste0("DownDar", 1:length(downs))
downsbedname <- paste0("downDARS_", experiment, ".bed")
export(downs, downsbedname)

# CHIPSEEKER ANNOTATIONS #######################

ups <- import.bed(upsbedname)
downs <- import.bed(downsbedname)

promoter <-
  getPromoters(TxDb = TxDb.Athaliana.BioMart.plantsmart28,
               upstream = 1000,
               downstream = 400)

tagMatrix_u <- getTagMatrix(ups, windows = promoter)
tagMatrix_d <- getTagMatrix(downs, windows = promoter)


pdf(file = paste0(experiment,"_TSS_heatmaps.pdf"), width = 6, height = 6)
par(mfrow = c(1, 2))
tagHeatmap(tagMatrix_u, xlim = c(-1000, 400), color = "green",title = "upDARs")
tagHeatmap(tagMatrix_d, xlim = c(-1000, 400), color = "red",title = "downDARs")
dev.off()

pdf(file = paste0(experiment,"_TSS_profiles.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))
plotAvgProf(
  tagMatrix_u,
  xlim = c(-1000, 400),
  conf = 0.95,
  resample = 1000,
  xlab = "Genomic Region (5'->3')",
  ylab = "Read Count Frequency Up DARs"
)

plotAvgProf(
  tagMatrix_d,
  xlim = c(-1000, 400),
  conf = 0.95,
  resample = 1000,
  xlab = "Genomic Region (5'->3')",
  ylab = "Read Count Frequency Down DARs"
)
dev.off()

peakAnno_u <- annotatePeak(ups,
                           tssRegion = c(-500, 200),
                           TxDb = TxDb.Athaliana.BioMart.plantsmart28,
                           annoDb = "org.At.tair.db")

peakAnno_d <- annotatePeak(
  downs,
  tssRegion = c(-500, 200),
  TxDb = TxDb.Athaliana.BioMart.plantsmart28,
  annoDb = "org.At.tair.db"
)

pdf(file = paste0(experiment,"_upDARs_pie.pdf"), width = 6, height = 6)
par(mfrow = c(1, 1))
plotAnnoPie(peakAnno_u)
dev.off()

pdf(file = paste0(experiment,"_downDARS_pie.pdf"), width = 6, height = 6)
par(mfrow = c(1, 1))
plotAnnoPie(peakAnno_d)
dev.off()


u <- as.data.frame(peakAnno_u)
d <- as.data.frame(peakAnno_d)

# Cluster files

b <-
  cbind(
    as.data.frame(merged_results$regions),
    as.data.frame(merged_results$best),
    as.data.frame(merged_results$combined)
  )

colnames(b) <- c(
  "seqnames",
  "start",
  "end",
  "width",
  "strand",
  "bestWindow",
  "best_logFC",
  "best_logCPM",
  "best_F",
  "best_PValue",
  "best_FDR",
  "nWindows",
  "logFC.up",
  "logFC.down",
  "combined_PValue",
  "combined_FDR",
  "direction"
)

write.table(
  b,
  file = paste0("allClusters_", experiment, ".tsv"),
  row.names = F,
  quote = F,
  sep = "\t"
)

# full tables

columnNames <- c(
  "seqnames",
  "start",
  "end",
  "width",
  "strand",
  "name",
  "score",
  "annotation",
  "geneChr",
  "geneStart",
  "geneEnd",
  "geneLength",
  "geneStrand",
  "geneID",
  "transcriptID",
  "distanceToTSS",
  "bestWindow",
  "best_logFC",
  "best_logCPM",
  "best_F",
  "best_PValue",
  "best_FDR",
  "nWindows",
  "logFC.up",
  "logFC.down",
  "combined_PValue",
  "combined_FDR",
  "direction"
)

b <-
  cbind(u,
        as.data.frame(temp1ups$best),
        as.data.frame(temp1ups$combined))
colnames(b) <- columnNames
write.table(
  b,
  file = paste0("upDARs_bigtable_", experiment, ".tsv"),
  row.names = F,
  quote = F,
  sep = "\t"
)

b <-
  cbind(d,
        as.data.frame(temp1dws$best),
        as.data.frame(temp1dws$combined))
colnames(b) <- columnNames
write.table(
  b,
  file = paste0("downDARs_bigatable_", experiment, ".tsv"),
  row.names = F,
  quote = F,
  sep = "\t"
)

# filtered results plot

pdf(file = paste0(experiment,"_filteredResults.pdf"), width = 6, height = 6) 

par(mfrow = c(1, 2))
plot(merged_results$best$logFC, main=paste("Unfiltered results"))
abline(h = 0.8 , col = "red", lwd = 2)
abline(h = -0.8 , col = "red", lwd = 2)

plot(tf1$best$logFC, main=paste("Filtered results"))
abline(h = 0.8 , col = "red", lwd = 2)
abline(h = -0.8 , col = "red", lwd = 2)

dev.off()

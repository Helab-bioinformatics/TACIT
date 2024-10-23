
######Creat gene activity matrix######

#remotes::install_version("Seurat", version = "3.1.2")
library(Seurat)
library(ggplot2)
library(patchwork)
library(pbapply)

CreateGeneActivityMatrix <- function (peak.matrix, annotation.file, seq.levels = chromosomes, include.body = TRUE, upstream = 2000, downstream = 0, 
                                      verbose = TRUE) 
{
  if (!PackageCheck("GenomicRanges", error = FALSE)) {
    stop("Please install GenomicRanges from Bioconductor.")
  }
  if (!PackageCheck("rtracklayer", error = FALSE)) {
    stop("Please install rtracklayer from Bioconductor.")
  }
  peak.df <- rownames(x = peak.matrix)
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, 
                                                            pattern = ":", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", "start", "end")
  #peak.df <- peak.df[order(peak.df$start),]
  #peak.df <- peak.df[-55611:-55629,]
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
  
  
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 
                                 0, ]) <- 1
  gtf <- rtracklayer::import(con = annotation.file)
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, 
                                     pruning.mode = "coarse")
  if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
    GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
  }
  gtf.genes <- gtf[gtf$type == "gene"]
  if (include.body) {
    gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, 
                            downstream = downstream)
  }
  else {
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, 
                                                     upstream = upstream, downstream = downstream)
  }
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, 
                                                     subject = gtf.body_prom)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance ==0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
  gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]
  #chose own gene type such as gene name（symbol） or gene id 
  peak.ids$gene.name <- gene.ids$gene_id
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
  annotations <- peak.ids[, c("peak", "gene.name")]
  colnames(x = annotations) <- c("feature", "new_feature")
  peak.matrix <- as(object = peak.matrix, Class = "matrix")
  all.features <- unique(x = annotations$new_feature)
  library(pbapply)
  mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
  
  newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x) {
    features.use <- annotations[annotations$new_feature == 
                                  all.features[[x]], ]$feature
    submat <- peak.matrix[features.use, ]
    if (length(x = features.use) > 1) {
      return(Matrix::colSums(x = submat))
    }
    else {
      return(submat)
    }
  })
  newmat <- t(x = newmat)
  rownames(x = newmat) <- all.features
  colnames(x = newmat) <- colnames(x = peak.matrix)
  return(as(object = newmat, Class = "dgCMatrix"))
}


Extend <- function (x, upstream = 0, downstream = 0) 
{
  if (any(GenomicRanges::strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- GenomicRanges::strand(x = x) == "+" | GenomicRanges::strand(x = x) == 
    "*"
  new_start <- GenomicRanges::start(x = x) - ifelse(test = on_plus, 
                                                    yes = upstream, no = downstream)
  new_end <- GenomicRanges::end(x = x) + ifelse(test = on_plus, 
                                                yes = downstream, no = upstream)
  IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, 
                                             end = new_end)
  x <- GenomicRanges::trim(x = x)
  return(x)
}

# setwd("/media/helab/data1/min/02_tacit/03_early_stage/early2cell/WNN/")
# 
cotacit_count <- read.csv("./cotacit.peaks.countMatrix.txt", sep = " ")
#remove non numeric rows
cotacit_count <- cotacit_count[grep("chrUn", rownames(cotacit_count), invert = T),]
cotacit_count <- cotacit_count[grep("random", rownames(cotacit_count), invert = T),]

annotation.file <- "/media/helab/data1/min/00_reference/Mus_musculus.GRCm39.110.gtf"
chromosomes = c(1:19, "X", "Y")
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = cotacit_count, annotation.file, 
                                            seq.levels = chromosomes, upstream = 5000,downstream = 5000,  verbose = TRUE)
geneActivity <- as.matrix(activity.matrix)


#######################################  transfer the rowname of activity.matrix to gene Symbol ###################################################################################

library(clusterProfiler)
library(tidyr)
library(dplyr)
library(org.Mm.eg.db)
library(openxlsx)
library(org.Mm.eg.db)

# eg2symbol=toTable(org.Mm.egSYMBOL)
# eg2name=toTable(org.Mm.egGENENAME)
# eg2ensembl=toTable(org.Mm.egENSEMBL)
# eg2alias=toTable(org.Mm.egALIAS2EG)

geneSymbol <-  bitr(rownames(geneActivity), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
head(geneSymbol)
dim(geneSymbol)
geneSymbol <- na.omit(geneSymbol)
length(unique(geneSymbol$SYMBOL))
geneSymbol <- geneSymbol[!duplicated(geneSymbol$SYMBOL),]
geneSymbol <- geneSymbol[!duplicated(geneSymbol$ENSEMBL),]
rownames(geneSymbol) <- geneSymbol$ENSEMBL
geneSymbol <- geneSymbol[rownames(geneActivity),]
rownames(geneActivity) <- geneSymbol$SYMBOL
geneActivity <- as.data.frame(geneActivity)
geneActivity <- geneActivity[!grepl("NA", rownames(geneActivity)),]

write.table(geneActivity, "./cotacit_tss5kb_ActivityMatrix.txt", quote = F)




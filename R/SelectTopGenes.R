#'  Select top or bottom N genes based on a selection criterion
#'
#' @description
#' The SelectTopGenes function enables selecting top or bottom N genes based on a criterion (e.g. log2FC or adj.p.value).
#'
#' @param gene.markers a data frame of the gene markers found by FindAllGeneMarkers function
#' @param top.N how many top or bottom genes to select (default 10)
#' @param criterion.type which criterion to use for selecting the genes (default "log2FC")
#' @param bottom whether to select bottom instead of top N genes (default FALSE)
#' @return data frame
#' @keywords select top bottom N genes
#' @importFrom dplyr group_by
#' @importFrom dplyr %>%
#' @importFrom dplyr top_n

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=15)
#' gene_markers <- FindAllGeneMarkers(iloreg_object,
#'                                   clustering.type = "manual",
#'                                   test = "wilcox",
#'                                   logfc.threshold = 0.25,
#'                                   min.pct = 0.1,
#'                                  min.diff.pct = NULL,
#'                                  pseudocount.use = 1,
#'                                  min.cells.group = 3,
#'                                 return.thresh = 0.01,
#'                                 only.pos = FALSE,max.cells.per.cluster = NULL)
#' ## Select top 10 and top 1 genes based on log2 fold-change and Bonferroni adjusted p-value.
#'top10_log2FC <- SelectTopGenes(gene_markers,top.N = 10,criterion.type = "log2FC",bottom = FALSE)
#'top1_log2FC <- SelectTopGenes(gene_markers,top.N = 1,criterion.type = "log2FC",bottom = FALSE)
#'top10_adj.p.value <- SelectTopGenes(gene_markers,top.N = 10,criterion.type = "adj.p.value",bottom = TRUE)
#'top1_adj.p.value <- SelectTopGenes(gene_markers,top.N = 1,criterion.type = "adj.p.value",bottom = TRUE)


SelectTopGenes <- function(gene.markers=NULL,top.N=10,criterion.type="log2FC",bottom=FALSE)
{

  if (bottom)
  {
    top.N <- -top.N
  }

  gene.markers %>% group_by(cluster) %>% top_n(top.N, get(criterion.type)) -> top_N

  return(top_N)
}

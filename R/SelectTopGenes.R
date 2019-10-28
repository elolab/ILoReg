#'  Select top or bottom N genes based on a selection criterion
#'
#' @description
#' The GeneDropoutRatePlot function enables visualizion of dropout rates (fraction of cells expressing a gene) against
#' the average non-zero expression values. This function can aid the user to determinine if a gene is differentially expressed or not,
#' since differentially expressed genes are likely to be above the expected curve. Cells can be filtered based on clustering before
#' calculating the values.
#'
#' @param gene.markers a character vector of the genes to visualize over the dropout rate curve
#' @param top.N return.plot whether to return the ggplot2 object or just draw it (default \code{FALSE})
#' @param criterion.type use data from these clusters only
#' @param bottom "manual" or "optimal". "manual" refers to the clustering formed using the "SelectKClusters" function and "optimal" to the clustering using the "CalculateSilhouetteInformation" function. (default "manual")
#' @return data frame
#' @keywords dropout rate curve non-zero average expression gene
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
#' GeneDropoutRatePlot(iloreg_object,genes = c("CD79A","TCL1A","IGLL5","VPREB3"),use.clusters = c(24,15,14,2),clustering.type = "manual")


SelectTopGenes <- function(gene.markers=NULL,top.N=10,criterion.type="log2FC",bottom=FALSE)
{

  if (bottom)
  {
    top.N <- -top.N
  }

  gene.markers %>% group_by(cluster) %>% top_n(top.N, get(criterion.type)) -> top_N

  return(top_N)
}

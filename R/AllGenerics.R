#' @export
setGeneric("PrepareILoReg",signature = "object",
           function(object) {
             standardGeneric("PrepareILoReg")
           })

#' @export
setGeneric("RunParallelICP",signature = "object",
           function(object, k = 15, d = 0.3, L = 200, r = 5, C = 0.3, reg.type = "L1", max.iter = 200, threads = 0,icp.batch.size=Inf) {
             standardGeneric("RunParallelICP")
           })

#' @export
setGeneric("RunPCA", signature = "object",
           function(object,
                    p = 50,
                    scale = FALSE,
                    threshold = 0) {
             standardGeneric("RunPCA")
           })

#' @export
setGeneric("PCAElbowPlot", signature = "object",
           function(object,
                    return.plot=FALSE) {
             standardGeneric("PCAElbowPlot")
           })

#' @export
setGeneric("RunUMAP", signature = "object",
           function(object) {
             standardGeneric("RunUMAP")
           })

#' @export
setGeneric("RunTSNE", signature = "object",
           function(object,
                    perplexity = 30) {
             standardGeneric("RunTSNE")
           })

#' @export
setGeneric("HierarchicalClustering", signature = "object",
           function(object) {
             standardGeneric("HierarchicalClustering")
           })

#' @export
setGeneric("CalcSilhInfo", signature = "object",
           function(object,
                    K.start = 2, K.end = 50) {
             standardGeneric("CalcSilhInfo")
           })

#' @export
setGeneric("SilhouetteCurve", signature = "object",
           function(object,
                    return.plot = FALSE) {
             standardGeneric("SilhouetteCurve")
           })

#' @export
setGeneric("SelectKClusters", signature = "object",
           function(object,
                    K = NULL) {
             standardGeneric("SelectKClusters")
           })

#' @export
setGeneric("MergeClusters", signature = "object",
           function(object,
                    clusters.to.merge = "",
                    new.name = "") {
             standardGeneric("MergeClusters")
           })

#' @export
setGeneric("RenameAllClusters", signature = "object",
           function(object,
                    new.cluster.names = "") {
             standardGeneric("RenameAllClusters")
           })

#' @export
setGeneric("RenameCluster", signature = "object",
           function(object,
                    old.cluster.name = "",
                    new.cluster.name = "") {
             standardGeneric("RenameCluster")
           })

#' @export
setGeneric("GeneScatterPlot", signature = "object",
           function(object,
                    genes = "",
                    return.plot = FALSE,
                    dim.reduction.type = "tsne",
                    point.size = 0.7,
                    title = "",
                    plot.expressing.cells.last = FALSE,
                    nrow = NULL,
                    ncol = NULL) {
             standardGeneric("GeneScatterPlot")
           })

#' @export
setGeneric("ClusteringScatterPlot",signature = "object",
           function(object,
                    clustering.type = "manual",
                    return.plot = FALSE,
                    dim.reduction.type = "",
                    point.size = 0.7,
                    title = "",
                    show.legend = TRUE) {
             standardGeneric("ClusteringScatterPlot")
           })

#' @export
setGeneric("FindAllGeneMarkers", signature = "object",
           function(object,
                    clustering.type = "manual",
                    test = "wilcox",
                    log2fc.threshold = 0.25,
                    min.pct = 0.1,
                    min.diff.pct = NULL,
                    min.cells.group = 3,
                    max.cells.per.cluster = NULL,
                    return.thresh = 0.01,
                    only.pos = FALSE) {
             standardGeneric("FindAllGeneMarkers")
           })

#' @export
setGeneric("FindGeneMarkers", signature = "object",
           function(object,
                    clusters.1 = NULL,
                    clusters.2 = NULL,
                    clustering.type = "",
                    test = "wilcox",
                    logfc.threshold = 0.25,
                    min.pct = 0.1,
                    min.diff.pct = NULL,
                    min.cells.group = 3,
                    max.cells.per.cluster = NULL,
                    return.thresh = 0.01,
                    only.pos = FALSE) {
             standardGeneric("FindGeneMarkers")
           })

#' @export
setGeneric("VlnPlot", signature = "object",
           function(object,
                    clustering.type = "manual",
                    genes = NULL,
                    return.plot = FALSE,
                    rotate.x.axis.labels=FALSE) {
             standardGeneric("VlnPlot")
           })

#' @export
setGeneric("GeneHeatmap", signature = "object",
           function(object,
                    clustering.type = "manual",
                    gene.markers = NULL) {
             standardGeneric("GeneHeatmap")
           })

#' @export
setGeneric("AnnotationScatterPlot", signature = "object",
           function(object,
                    annotation = NULL,
                    return.plot = FALSE,
                    dim.reduction.type = "",
                    point.size = 0.7,
                    show.legend = FALSE) {
             standardGeneric("AnnotationScatterPlot")
           })


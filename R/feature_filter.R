
# Filter out features with low counts
# Uses raw sequence read counts and number of samples per group to filter 
# out genes with low number of read counts in most samples.
#' @param rawcounts counts matrix
#' @param n_samples minimum number of samples passing the filter
#' @return logical vector
#' @import edgeR
#' @export
feature_filter <- function(rawcounts, n_samples) {
  
  # Calculate normalised counts
  cpms <- edgeR::cpm(rawcounts)
  
  # Smallest libsize in millions of counts
  libsize <- min(colSums(rawcounts, na.rm = TRUE)) / 1e6
  
  # Library size adjusted feature filter
  rowSums(cpms > 10 / libsize) >= n_samples
}

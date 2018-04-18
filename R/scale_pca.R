#' Scale function for PCA model
#' 
#' @description Before the data is used for PCA, the data should be centered and scaled
#' @param normal The standard matrix is used to scale
#' @param fault The matrix is used to be scaled. W
#' @param ... Lazy dots for additional internal arguments
#' 
#' @return A matrix of the scaled matrix fault
#' 
#' @examples 
#' data(normal_data)
#' normal_scale <- scale_pca(normal_data)
#' @export
scale_pca <- function(normal, fault,...){
  normal_mean <- apply(normal,2,mean)
  normal_sd <- apply(normal,2,sd)
  
  fault_scale <- apply(fault, 1, function(x,normal_mean,normal_sd){
    (x-normal_mean)/normal_sd
  },normal_mean,normal_sd)
  return(t(fault_scale))
}
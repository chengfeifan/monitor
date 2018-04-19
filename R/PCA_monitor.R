#' PCA model for process monitor
#' 
#' @description Calculate the principal component analysis for process monitor,
#' and also find the squared prediction error (SPE) and Hotelling's T2 test
#' statistic values for each observation in this data matrix.
#' 
#' @param data A centered-and-scaled data matrix
#' @param kernel_num The number of principle component
#' @param ... Lazy dots for additional internal arguments
#'
#' @return A list of class "pca" with the following:
#' \itemize{
#'   \item{projectionMatrix -- }{the q eigenvectors corresponding to the q
#'     largest eigenvalues as a p x q projection matrix}
#'   \item{LambdaInv -- }{the diagonal matrix of inverse eigenvalues}
#'   \item{SPE -- }{the vector of SPE test statistic values for each of the n
#'     observations contained in "data"}
#'   \item{T2 -- }{the vector of Hotelling's T2 test statistic for each of the
#'     same n observations}
#' }
#' @references https://github.com/gabrielodom/mvMonitoring
#' @examples
#' 
#' @export
#' 
PCA_model <- function(data,kernel_num,...){
  UseMethod('PCA_model')
}

#' @keywords internel
#' @export
#' @importFrom stats cor
PCA_model.matrix <- function(data,kernel_num,...){
  R <- cor(data, use = "pairwise.complete.obs")
  eigenR <- eigen(R)
  evalR <- eigenR$values
  evecR <- eigenR$vectors
  P <- as.matrix(evecR[, 1:kernel_num]) # Transformation matrix
  Lambda <- diag(evalR[1:kernel_num], ncol = length(1:kernel_num))
  
  # Rotated Matrix
  PCs <- data %*% P
  
  # Reduced Matrix in Original Space
  X.hat <- PCs %*% t(P)
  
  # Residual Matrix
  E <- data - X.hat
  
  # Squared prediction error monitoring statistic
  SPEs <- diag(E %*% t(E))
  
  # Hotelling's T^2 monitoring statistic
  LambdaInv <- solve(Lambda)
  T2s <- diag(PCs %*% LambdaInv %*% t(PCs))
  
  object <- list(projectionMatrix = P,
                 LambdaInv = LambdaInv,
                 SPE = SPEs,
                 T2 = T2s)
  
  class(object) <- "PCA_model"
  return(object)
}


#' Calculate threhold for PCA model
#' @description Calculate the non-parametric critical value threshold estimates
#' for the SPE and T2 monitoring test statistics.
#' 
#' @param pca_object A list with class "PCA_model" from the internal PCA_model() function
#' @param alpha The upper 1 - alpha quantile of the SPE and T2 densities from
#'   the training data passed to this function. Defaults to 0.01.
#' @param ... Lazy dots for additional internal arguments
#' @return A list with classes "threshold" and "pca" containing:
#'   \itemize{
#'     \item{SPE_threshold -- }{the 1 - alpha quantile of the estimated SPE
#'       density}
#'     \item{T2_threshold -- }{the 1 - alpha quantile of the estimated Hotelling's
#'       T2 density}
#'     \item{projectionMatrix -- }{a projection matrix from the data feature space
#'       to the feature subspace which preserves some pre-specified proportion
#'       of the energy of the data scatter matrix.}
#'     \item{LambdaInv -- }{a diagonal matrix of the reciprocal eigenvalues of the
#'       data scatter matrix}
#'     \item{T2 -- }{the vector of Hotelling's T2 test statistic values for each of the n
#'     observations in "data"}
#'     \item{SPE -- }{the vector of SPE test statistic values for each of the n
#'     observations in "data"}
#' 
#'   }
#' @details This function takes in a pca object returned by the pca() function
#' and a threshold level defaulting to alpha = 0.1 percent of the
#' observations. This critical quantile is set this low to reduce false
#' alarms, as described in Kazor et al (2016). The function then returns a
#' calculated SPE threshold corresponding to the 1 - alpha critical value, a
#' similar T2 threshold, and the projection and Lambda Inverse (1 /
#' eigenvalues) matrices passed through from the pca() function call.
#' @export
threshold <- function(pca_object, alpha = 0.01, ...){
  UseMethod("threshold")
}

#' @export
#' @keywords internal
#'
#' @importFrom BMS quantile.density
#' @importFrom stats density
#'
threshold.PCA_model <- function(pca_object, alpha = 0.01, ...){
  
  spe <- pca_object$SPE
  t2 <- pca_object$T2
  
  SPE.np.dens <- density(spe,
                         bw = "SJ", # Sheather Jones
                         kernel = "gaussian",
                         from = 0)
  
  # BMS::quantile.density
  SPE.lim.np <- quantile.density(SPE.np.dens, 1 - alpha)
  
  T2.np.dens <- density(t2,
                        bw = "SJ", # Sheather Jones
                        kernel = "gaussian",
                        from = 0)
  # BMS::quantile.density
  T2.lim.np <- quantile.density(T2.np.dens, 1 - alpha)
  
  object <- list(SPE_threshold = SPE.lim.np,
                 T2_threshold = T2.lim.np,
                 projectionMatrix = pca_object$projectionMatrix,
                 LambdaInv = pca_object$LambdaInv,
                 T2 = pca_object[[4]],
                 SPE = pca_object[[3]])
  
  class(object) <- c("threshold", "PCA_model")
  object
}

#' Process detection 
#' @description After the model builds, use it to detect if a single
#'  multivariate observation is beyond normal operating conditions.
#'  
#' @param threshold_object An object of classes "threshold" and "pca" returned
#'   by the internal threshold() function.
#' @param observation A single row of an xts data matrix  to
#'   compare against the thresholds
#' @param ... Lazy dots for additional internal arguments
#'
#' @return A named 1 x 4 matrix with the following entries for the single row
#'   observation passed to this function:
#'   \itemize{
#'     \item{SPE -- }{the SPE statistic value}
#'     \item{SPE_Flag -- }{the SPE fault indicator, where 1 represents a flag and
#'       0 marks that the observation is within the normal operating conditions}
#'     \item{T2 -- }{the T2 statistic value}
#'     \item{T2_Flag -- }{the T2 fault indicator, defined the same as SPE_Flag}
#'   }
#'
#' @details This function takes in a threshold object returned by the
#'   threshold() function and a single observation as a matrix or xts row.
#'   Internally, the function multiplies the observation by the projection
#'   matrix returned within the threshold object, calculates the SPE and T2
#'   process monitoring statistics for that observation, and compares these
#'   statistics against their corresponding threshold values to determine if the
#'   observation lies outside the expected boundaries. The function then returns
#'   a row vector of the SPE test statistic, a logical indicator marking if this
#'   statistic is beyond the threshold, the Hotelling's T2 statistic, and an
#'   indicator if this statistic is beyond the threshold. Observations with
#'   monitoring statistics beyond the calculated threshold are marked with a 1,
#'   while observations within normal operating conditions are marked with a 0.
#'   These threshold values are passed from the threshold() function through
#'   this function via a returned threshold object. This object will be used in
#'   higher function calls.
#' @examples 
#' data(normal_data)
#' data(fault_data)
#' normal_embed <- embed(scale_pca(normal = normal_data, fault = normal_data),2)
#' fault_embed <- embed(scale_pca(normal=normal_data,fault=fault_data),2)
#' model_pca <- PCA_model(data=normal_embed, kernel_num=5)
#' threshold_pca <- threshold(pca_object=model_pca)
#' fault_detection <- faultDetect(threshold_object = threshold_pca, observation =  fault_embed)
#' @export
faultDetect <- function(threshold_object, observation, ...){
  UseMethod("faultDetect")
}

#' @export
#' @keywords internal
#'
faultDetect.threshold <- function(threshold_object, observation, ...){
  
  SPEthreshold <- threshold_object$SPE_threshold
  T2threshold <- threshold_object$T2_threshold
  P <- threshold_object$projectionMatrix
  LambdaInv <- threshold_object$LambdaInv
  
  proj_observation <- observation %*% P
  
  # Reduced Observation in Original Space
  obs.hat <- proj_observation %*% t(P)
  
  # Residual Vector
  E <- observation - obs.hat
  
  # Squared prediction error monitoring statistic
  SPE <- diag(E %*% t(E))
  SPE_flag <- as.numeric(SPE > SPEthreshold)
  
  # Hotelling's T^2 monitoring statistic
  T2 <- diag(proj_observation %*% LambdaInv %*% t(proj_observation))
  T2_flag <- as.numeric(T2 > T2threshold)
  object <- list(SPE=SPE,
                           SPE_Flag=SPE_flag,
                           SPE_threshold=SPEthreshold,
                           T2=T2,
                           T2_Flag=T2_flag,
                           T2_threshold=T2threshold)
  class(object)<-c('threshold','faultDetect')
  object
}

#' Plot for the class threshold
#' @description Plot the statistics T2 and SPE of class threshold.
#' 
#' @param threshold_object The class of threshold
#' 
#' @return A plot
#' @export
#' @importFrom gridExtra grid.arrange
plot.threshold <- function(threshold_object){
  p1 <- qplot(x=c(1:length(threshold_object$SPE)),threshold_object$SPE,geom = 'line') + theme_classic() + xlab('Time') + ylab('SPE') +
    geom_hline(yintercept = threshold_object$SPE_threshold,col='red',lty=2)
  p2 <- qplot(x=c(1:length(threshold_object$T2)),threshold_object$T2,geom = 'line') + theme_classic() + xlab('Time') + ylab('T2') +
    geom_hline(yintercept = threshold_object$T2_threshold,col='red',lty=2)
  grid.arrange(p1, p2, ncol = 1)
}
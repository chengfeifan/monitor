#' kernel PCA for process monitor
#' 
#' @description Use kernel principle component analysis for process monitoring
#' and also find the squared prediction error (SPE) and Hotelling's T2 test
#' statistic values for each observation in this data matrix.
#' 
#' @param data A centered-and-scaled data matrix
#' @param kernel The kernel function to be used to calculate the kernel matrix. 
#' This has to be a function of class kernel,
#'i.e. which can be generated either one of the build in kernel generating 
#'functions (e.g., rbfdot etc.) or a user defined function of 
#'class kernel taking two vector arguments and returning a scalar.
#'
#' @param kernel_num The number of principle component
#' @param ... Lazy dots for additional internal arguments
#' @references http://modeleau.fsg.ulaval.ca/fileadmin/modeleau/documents/Publications/pvr487.pdf
#' @return A list of class 'KPCA' with following
#' \itemize{
#'   \item{projectionMatrix -- }{the q eigenvectors corresponding to the q
#'     largest eigenvalues as a p x q projection matrix}
#'   \item{LambdaInv -- }{the diagonal matrix of inverse eigenvalues}
#'   \item{SPE -- }{the vector of SPE test statistic values for each of the n
#'     observations contained in "data"}
#'   \item{T2 -- }{the vector of Hotelling's T2 test statistic for each of the
#'     same n observations}
#'   \item{model_data -- }{the data used to build the model}
#'   \item{K_hat -- }{The mean and center matrix of kernel matrix}
#'   \item{eigenK -- }{A list consists of eigen values and eigen vector of K_hat}
#' }
#' @example
#' @export
KPCA_model <- function(data,kernel,kernel_num,...){
  UseMethod('KPCA_model')
}

#' @keywords 
#' @export
#' @importFrom kernlab kernelMatrix
KPCA_model <- function(data,kernel,kernel_num,...){
  K <- kernelMatrix(kernel = kernel,data)@.Data
  m <- ncol(K)
  # center and mean
  # K_hat = K - 1^m*K - K*1^m + 1^m *K* 1^m
  # 1^m is a matrix of m*m, whose element is 1/M
  M = 1/m * matrix(1,ncol=m,nrow = m)
  K_hat = K - M %*% K - K %*% M + M %*% K %*% M
  
  # eigen value of K_hat
  eigenK <- eigen(K_hat)
  evalK <- eigenK$values
  evecK <- eigenK$vectors
  
  P <- Re(as.matrix(evecK[, 1:kernel_num]))
  
  # Transformation matrix
  Lambda <- Re(diag(evalK[1:kernel_num], ncol = length(1:kernel_num)))
  
  # Rotated Matrix
  PCs <- K_hat %*% P
  
  # Hotelling's T^2 monitoring statistic
  LambdaInv <- solve(Lambda)
  T2s <- diag(PCs %*% LambdaInv %*% t(PCs))
  
  # Calculate Squared prediction error monitoring statistic
  PCs_all <- Re(K_hat %*% evecK)
  SPEs <- diag(PCs_all %*% t(PCs_all)) -  diag(PCs %*% t(PCs))
  
  object <- list(projectionMatrix = P,
                 LambdaInv = LambdaInv,
                 SPE = SPEs,
                 T2 = T2s,
                 model_data=data,
                 K=K,
                 eigenK=eigenK,
                 kernel_num=kernel_num,
                 kernel = kernel)
  
  class(object) <- "KPCA_model"
  return(object)
}

#' @export
#' @keywords internal
#'
#' @importFrom BMS quantile.density
#' @importFrom stats density
#'
threshold.KPCA_model <- function(kpca_object, alpha = 0.01,empirical = TRUE, ...){
  spe <- kpca_object$SPE
  t2 <- kpca_object$T2
  
  if(empirical){
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
  } else{
    p = kpca_object$kernel_num
    N = nrow(kpca_object$model_data)
    
    # T2 limit
    T2.lim.np <- p*(N-1)/(N-p)*qf(p=1-alpha,df1 = p,df2 = N)
    
    # SPE limit
    a <- mean(spe)
    b <- var(spe)
    g <- b/(2*a)
    h <- 2*a^2/b
    SPE.lim.np <- g*qchisq(p=1-alpha,df=h)
  }
  
  object <- list(SPE_threshold = SPE.lim.np,
                 T2_threshold = T2.lim.np,
                 projectionMatrix = kpca_object$projectionMatrix,
                 LambdaInv = kpca_object$LambdaInv,
                 T2 = kpca_object$T2,
                 SPE = kpca_object$SPE,
                 model_data = kpca_object$model_data,
                 K=kpca_object$K,
                 eigenK=kpca_object$eigenK,
                 kernel_num=kpca_object$kernel_num,
                 kernel = kpca_object$kernel)
  
  class(object) <- c("KPCA_model", "threshold")
  object
}

#' @export
#' @keywords internal
#' @importFrom kernlab kernelMatrix
faultDetect.KPCA_model <- function(threshold_object, observation, ...){
  kernel <- threshold_object$kernel
  model_data <- threshold_object$model_data
  SPEthreshold <- threshold_object$SPE_threshold
  T2threshold <- threshold_object$T2_threshold
  
  L = nrow(observation)
  m = nrow(model_data)
  
  M=1/m*matrix(1,nrow=m,ncol=m)
  t_M = 1/m*matrix(1,nrow=L,ncol=m)
  
  K_test <- kernelMatrix(kernel = kernel,x=observation,y=model_data)@.Data
  K <- threshold_object$K
  # K_test = K_test - 1'M * K - K_test * 1M + 1'M * K * 1M
  # 1'M is L*K matrix, whose element is 1/M
  K_test_hat <- K_test -  t_M%*% K - K_test %*% M + t_M %*% K %*% M
  
  PCs <- Re(K_test_hat %*% threshold_object$projectionMatrix)
 
  # Hotelling's T^2 monitoring statistic
  T2 <- diag(PCs %*% threshold_object$LambdaInv %*% t(PCs))  
  T2_flag <- as.numeric(T2 > T2threshold)
  
  # Calculate Squared prediction error monitoring statistic
  PCs_all <- Re(K_test_hat %*% threshold_object$eigenK$vectors)
  SPE  <- diag(PCs_all %*% t(PCs_all)) -  diag(PCs %*% t(PCs))
  SPE_flag <- as.numeric(SPE > SPEthreshold)
  
  object <- list(SPE=SPE,
                 SPE_Flag=SPE_flag,
                 SPE_threshold=SPEthreshold,
                 T2=T2,
                 T2_Flag=T2_flag,
                 T2_threshold=T2threshold)
  class(object)<-c('threshold','faultDetect')
  object
}

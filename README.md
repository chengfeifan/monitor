# monitor
Write the algorithm for process monitoring. The method is *PCA*, *KPCA*, *DPCA* and *DKPCA*. 

## PCA
Run PCA on a linear model

  library(monitor)
  data(normal_data)
  data(fault_data)
  normal_embed <- scale_pca(normal = normal_data, fault = normal_data)
  fault_embed <- scale_pca(normal=normal_data,fault=fault_data)
  model_pca <- PCA_model(data=normal_embed, kernel_num=2)
  threshold_pca <- threshold(pca_object=model_pca)
  fault_detection <- faultDetect(threshold_object = threshold_pca, observation =  fault_embed)

The Result figure is as following

![PCA monitoring]('./figure/PCA.png')

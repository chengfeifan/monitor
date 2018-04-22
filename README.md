# monitor
Write the algorithm for process monitoring. The method is *PCA*, *KPCA*, *DPCA* and *DKPCA*. 

## PCA
Run PCA on a linear autocorrelated model
```
library(monitor)
data(normal_data)
data(fault_data)
normal_embed <- scale_pca(normal = normal_data, fault = normal_data)
fault_embed <- scale_pca(normal=normal_data,fault=fault_data)
model_pca <- PCA_model(data=normal_embed, kernel_num=2)
threshold_pca <- threshold(pca_object=model_pca)
fault_detection <- faultDetect(threshold_object = threshold_pca, observation =  fault_embed)
```
The figure shows the result as follows:

![PCA monitoring](/figure/PCA.png)

## DPCA
Run DPCA on a linear autocorrelated model
```
# DPCA
library(monitor)
data(normal_data)
data(fault_data)
normal_embed <- embed(scale_pca(normal = normal_data, fault = normal_data),2)
fault_embed <- embed(scale_pca(normal=normal_data,fault=fault_data),2)
model_pca <- PCA_model(data=normal_embed, kernel_num=5)
threshold_pca <- threshold(pca_object=model_pca)
fault_detection <- faultDetect(threshold_object = threshold_pca, observation =  fault_embed)
plot(fault_detection)
```
The figure shows the result as follows:

![DPCA monitoring](/figure/DPCA.png)

## KPCA
Run KPCA on a nonlinear model

```
data("nonlinear")
nonlinear_train <- nonlinear[1:400,]
nonlinear_test <- nonlinear[401:nrow(nonlinear),]
normal_embed <- scale_pca(normal = nonlinear_train, fault = nonlinear_train)
fault_embed <- scale_pca(normal=nonlinear_train,fault=nonlinear_test)

kernel = rbfdot(sigma=1/7200)
model_kpca <- KPCA_model(data = normal_embed, kernel = kernel, kernel_num = 2)
threshold_kpca <- threshold(model_kpca)
fault_detection <- faultDetect(threshold_object = threshold_kpca, observation =  fault_embed)
plot(fault_detection)
```
The figure shows the result as follows:

![KPCA monitoring](/figure/KPCA.png)

## DKPCA
Run DKPCA on a nonlinear model

```
data("nonlinear")
nonlinear_train <- nonlinear[1:400,]
nonlinear_test <- nonlinear[401:nrow(nonlinear),]
normal_embed <- embed(scale_pca(normal = nonlinear_train, fault = nonlinear_train),2)
fault_embed <- embed(scale_pca(normal=nonlinear_train,fault=nonlinear_test),2)

kernel = rbfdot(sigma=1/7200)
model_kpca <- KPCA_model(data = normal_embed, kernel = kernel, kernel_num = 5)
threshold_kpca <- threshold(model_kpca)
fault_detection <- faultDetect(threshold_object = threshold_kpca, observation =  fault_embed)
plot(fault_detection)
```
The figure shows the result as follows:

![DKPCA monitoring](/figure/DKPCA.png)

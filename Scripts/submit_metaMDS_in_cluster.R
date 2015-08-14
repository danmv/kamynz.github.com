#!/usr/local/bin/Rscript --save 

## NOTES:

## Principal Coordinate Analysis => METRIC MULTIDIMENSIONAL SCALING OR classical scaling
## Objective:  transforms a distance matrix into a set of coordinates 
## such that the (Euclidean) distances derived from these coordinates approximate as well as possible the original distances

# install.packages(c("ecodist", "labdsv", "ape", "ade4", "smacof"))

#in order to use these packages => the matrix have to be converted to dist class

## loading en .rdata from cmdscale before
load("/biologia-scratch3/mc.martinez297/TESIS/semana1/Doing_Non_Metric_Multidimensional_Scaling/data_for_MDS.rdata")

# matrix to use the packages
matrix_for_MDS_Ecoli<-matrix_for_PCA_Ecoli

matrix_for_MDS_normalized_by_max_Ecoli<-matriz_Ecoli_for_PCA_normalized

matrix_for_MDS_normalized_by_colsums_Ecoli<-matriz_for_PCA_normalized_colSums_Ecoli

#View(matrix_for_MDS)

## getting distance matrixes in order to use the functions 

#dist matrix for MDS from non-normalized input data Ecoli
dist_matrix_for_MDS_Ecoli<-as.dist(matrix_for_MDS_Ecoli)

#dist matrix for MDS from normalized by max input data Ecoli
dist_for_MDS_normalized_by_max_Ecoli<-as.dist(matrix_for_MDS_normalized_by_max_Ecoli)

#dist matrix for MDS from normalized by colsums input data Ecoli
dist_for_MDS_normalized_by_colsums_Ecoli<-as.dist(matrix_for_MDS_normalized_by_colsums_Ecoli)



# loading all libraries 
library(vegan)
# library(ecodist)
# library(labdsv)
# library(ape)
# library(ade4)
# library(smacof)

## using the different function and seeing how much it takes for each one

## 1. Checking distance matrix

## 2. Using cmdscale() to perform a classical scaling 

## DOING CMDSCALE CLASICO

# # se para cuando le pongo mis datos
# system.time(MDS_from_cmdscale_non_normalized_Ecoli_input<-cmdscale(dist_matrix_for_MDS,k=3))
# 

#&&&&&&&&&&&&&&&&&&&&&&&& 1 &&&&&&&&&&&&&&&&&&&&&&&&#

#Using metaMDS from vegan to do Non-metric multidimensional scaling 

# data(dune)

# para saber si son iguales

# all(matrix_for_MDS_Ecoli==matrix_for_MDS_normalized_by_colsums_Ecoli)
# 
# all(matrix_for_MDS_normalized_by_colsums_Ecoli==matrix_for_MDS_normalized_by_max_Ecoli)

## IMPORTANTE

"The metaMDS function was designed to be used with community data. 
If you have other type of data, you should probably set some arguments 
to non-default values: probably at least wascores, autotransform and 
noshare should be FALSE. If you have negative data entries, metaMDS will 
set the previous to FALSE with a warning"

### Doing NMDS of Ecoli non-normalized

#min(dist_matrix_for_MDS_Ecoli)


Output_metaMDS_Ecoli_non_normalized<-metaMDS(comm = dist_matrix_for_MDS_Ecoli,distance="bray",
                                             k=3,wascores=F,autotransform=F,noshare=F)


Output_metaMDS_Ecoli_normalized_by_max<-metaMDS(comm = dist_for_MDS_normalized_by_max_Ecoli,distance="bray",
                                                k=3,wascores=F,autotransform=F,noshare=F)

Output_metaMDs_Ecoli_normalized_by_ColSums<-metaMDS(comm = dist_for_MDS_normalized_by_colsums_Ecoli,distance="bray",
                                                    k=3,wascores=F,autotransform=F,noshare=F)

## saving the results

results_MDS_from_Ecoli<-c("Output_metaMDS_Ecoli_non_normalized","Output_metaMDS_Ecoli_normalized_by_max",
                          "Output_metaMDs_Ecoli_normalized_by_ColSums")

## exporting the results in .rdata file

save(list=results_MDS_from_Ecoli,file="/biologia-scratch3/mc.martinez297/TESIS/semana1/Doing_Non_Metric_Multidimensional_Scaling/results_MDS_from_Ecoli.rdata")

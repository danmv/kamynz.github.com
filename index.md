---
layout: post
title: First result
description: ""
category: Thesis
tags: [thesis,first image]
---
{% include JB/setup %}

##Creating an image with PCA results

This is the first attempt of creating an image with dimension reduction using PCA. However, the sum of the variance 
of the first three components is ~20%. 

## R script with 5 functions to generate the image

    # GETTING LIBRARIES
    
    library("seqinr")
    library("ade4")
    library("Biostrings")
    
    
    ## funcion 1
    
    Processing_fragments<-function(PATH_FILES){
      
      #GETTING THE FILES AFTER FRAGMENTS OF 500
      files <- list.files(path=PATH_FILES, pattern=".fna500mer", full.names=T, recursive=FALSE)
      
      #GETTING THE FILES READING AS FASTA
      ncrna <- lapply(files, function(x) { read.fasta(x,seqonly = T) })
      
      
      fragmentsGeno1<-list()
      for(k in seq_along(ncrna[1]))
      {
        for(l in 1:10484)
        {
          fragmentsGeno1[l]<-ncrna[[k]][[l]]
          
        }
      }
      
      fragmentsGeno2<-list()
      for(k in seq_along(ncrna[2]))
      {
        for(l in 1:length(ncrna[[2]]))
        {
          fragmentsGeno2[l]<-ncrna[[k]][[l]]
          
        }
      }
    
      #GETTING ALL FRAGMENTS
      
      allFragments<-c(fragmentsGeno1,fragmentsGeno2)
    
      return(allFragments)
    
    }
    
    
    ## funcion 2
    
    Getting_frequency_account<-function(allFragments,kmer){
      
      #CONVERTING LOS FRAGMENTOS DE CADA FILE A OBJETOS DE DNAString
      
      DNA_String_Set_list_ALL<-list()
      
      for(i in seq_along(allFragments))
      {
        DNA_String_Set_list_ALL[i]<-DNAStringSet(allFragments[[i]])
      }
      
      # counting oligonucleotide
      countGenome1_Tetra<-lapply(DNA_String_Set_list_ALL,function(x) {oligonucleotideFrequency((x),kmer, as.prob = T) })
      
      # MATRIX FOR THE PCA
      
      #names columns
      col_names<-dimnames(countGenome1_Tetra[[1]])
      col_names<-col_names[[2]]
      
      #names rows
      frag_names<-c(paste("frag",c(1:length(allFragments)),sep=""))
      
      #matrix for PCA
      matrix_PCA<-matrix(unlist(countGenome1_Tetra),nrow = length(allFragments),ncol=256,byrow = T,dimnames=list(frag_names,col_names))
      
      return(matrix_PCA)
      
    }
    
    
    # View(matrix_PCA)
    
    
    ## funcion 3
    
    Getting_first_three_components<-function(matrix_PCA){
    
      ######## PCA with prcomp#########
      
      prcomp_All<-prcomp(matrix_PCA)
      
      #obtaing the sum of varianza of the first three components
      
      Var<-prcomp_All$sdev^2 / sum(prcomp_All$sdev^2)
      
      Varianza_3_first_comp<-Var[1:3]
      
      Varianza_3_first_comp_Porcent<-Varianza_3_first_comp*100
      
      Suma_total<-sum(Varianza_3_first_comp_Porcent)
      
      ## obteniendo eigen of first three components 
      
      loadings_prcomp<-prcomp_All$x
      
      #dim(loadings_prcomp)
      
      First_three_components<-loadings_prcomp[,c(1,2,3)]
      
      return(First_three_components)
    
    }
    
    #funcion 4
    
    Generating_hex_color_codes<-function(First_three_components){
      
      # getting min and max
      min<-min(First_three_components)
      max<-max(First_three_components)
      
      # getting ranges
      range_2_color<-c(min,max)
      range_RGB_color<-c(0,1)
      
      #making linear regression
      lm.out<-lm(range_RGB_color~range_2_color)
      
      #getting slope and intercept
      slope<-lm.out$coefficients[2]
      intercept<-lm.out$coefficients[1]
      
      #normalizing pca results to RGB
      new_Matriz<-(First_three_components*slope)+intercept
      
      new_Matriz<-as.matrix(new_Matriz)
      
      #using funcion rgb to generate matrix of hex color code
      
      #hex_Color_Matriz<-t(mapply(rgb, split(new_Matriz[,1], new_Matriz[,2],new_Matriz[,3],maxColorValue=255)))
      
      hex_Color_Vector<-vector()
      
      # list de cada r,g,b de cada fragmento
      
      rgb_List_Each_Fragment<-list()
      
      row_Final<-length(new_Matriz[,1])
      
      columns_Final<-length(new_Matriz[1,])
      
      for(i in 1:row_Final){
        
        for(j in 1:columns_Final){
          
          red<-new_Matriz[i,1]
          green<-new_Matriz[i,2]
          blue<-new_Matriz[i,3]
          
          hex_Color_Vector[i]<-rgb(red,green,blue,maxColorValue = 1)
          
          rgb_List_Each_Fragment[i]<-list(c(red,green,blue))
          
        }
        
      }
      
      return(rgb_List_Each_Fragment)
      
    }
    
    # Calling all the funcionts in order
    
    allFragments<-Processing_fragments("/Users/CamilaMV/Desktop/TESIS")
    
    matrix_PCA<-Getting_frequency_account(allFragments,4)
    
    First_three_components<-Getting_first_three_components(matrix_PCA)
    
    Hex_color_list<-Generating_hex_color_codes(First_three_components)
    
    
    head(Hex_color_list)
    
    
    ## CREANDO IMAGEN, SERA QUE TENGO QUE HACER ESCALA DE COLORES?
    
    ## test 3
    
    hex_10710<-hex[1:10710]
    
    length(hex_10710)
    
    hex_matriz_10710<-matrix(unlist(hex_10710), nrow=102, ncol = 105, byrow = TRUE)
    
    
    dim(hex_matriz_10700)
    
    #creating image
    
    library(grid)
    
    grid.raster(hex_matriz_10710, interpolate=FALSE)
    
    
    ## FUNCIONES ADICIONALES
    
    #finding factors of a number => brute force
    
    FUN <- function(x) {
      x <- as.integer(x)
      div <- seq_len(abs(x))
      factors <- div[x %% div == 0L]
      factors <- list(neg = -factors, pos = factors)
      return(factors)
    }
    
    FUN(10710)

## To-Do

- `I have to normalize` the input data for the count of oligonucleotide frequency of Escheria coli. 
    








# funciones que han pasado varias revisiones 

# GETTING LIBRARIES

library("seqinr")
library("ade4")
library("Biostrings")


## funcion 1


Processing_fragments<-function(PATH_FILES){
  
  #GETTING THE FILES AFTER FRAGMENTS OF 500
  #files <-list.files(path=PATH_FILES, pattern=".fna_2025_sesgos", full.names=T, recursive=FALSE)
  
  files <-list.files(path=PATH_FILES, pattern=".fna_2025_fragments", full.names=T, recursive=FALSE)
  
  #GETTING THE FILES READING AS FASTA
  ncrna<-lapply(files, function(x) { read.fasta(x,as.string = T,seqonly = T) })
  
  ##GETTING ALL THE FRAGMENTS
  
  fragments_length_list<-vector()
  
  ## Loop to get the number of fragments for each file
  for(i in 1:length(files))
  {
    n<-length(ncrna[[i]])
    fragments_length_list[i]<-as.numeric(n)
    
  }
  
  ## APPENDING EACH TEMP LIST OF FRAGMENTS IF ALLFRAGMENTS
  
  allFragmentsBetter<-list()
  
  for(j in 1:length(ncrna))
  {
    inicio = 1
    final = fragments_length_list[j]
    
    #cat(inicio,"\n",final,"\n")
    
    temp_genome_list<-list()
    
    for(k in 1:final)
    {
      temp_genome_list[k]<-ncrna[[j]][[k]]     
    }
    
    allFragmentsBetter<-append(allFragmentsBetter,temp_genome_list)
    
  }
  
  return(allFragmentsBetter)
  
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

# Getting_first_three_components_MODI<-function(matrix_for_PCA){
#   
#   ######## PCA with prcomp#########
#   
#   prcomp_All<-prcomp(matrix_for_PCA)
#   
#   #obtaing the sum of varianza of the first three components
#   
#   Var<-prcomp_All$sdev^2 / sum(prcomp_All$sdev^2)
#   
#   Variance_all_3<-Var[1:3]
#   
#   Suma_total<-sum((Variance_all_3)*100)
#   
#   ## obteniendo eigen of first three components 
#   
#   loadings_prcomp<-prcomp_All$x
#   
#   Main_3_Comps<-loadings_prcomp[,c(1,2,3)]
#   
#   #dim(loadings_prcomp)
#   
#   return(list("Var_C1_Red"=Variance_all_3[1],"Var_C2_Green"=Variance_all_3[2],
#               "Var_C3_Blue"=Variance_all_3[3],"Total_First_3_Comps"=Suma_total,
#               "Main_3_Comps"=Main_3_Comps,"pca"=prcomp_All))
#   
# }



#funcion 4

Generating_Hex_Vector_MODI<-function(Main_3_Comps){
  
  # getting min and max for the vector of each color
  
  minAll<-apply(Main_3_Comps,2,min)
  maxAll<-apply(Main_3_Comps,2,max)
  #   
  #   Main_3_Comps
  
  ## getting ranges => 1=Red, 2=Blue, 3=Green
  range_2_color_RED<-c(minAll[1],maxAll[1])
  
  range_2_color_GREEN<-c(minAll[2],maxAll[2])
  
  range_2_color_BLUE<-c(minAll[3],maxAll[3])
  
  range_RGB_color<-c(0,1)
  
  #   #making linear regression
  #   lm_all.out<-c(lm(range_RGB_color~range_2_color_RED),
  #                 lm(range_RGB_color~range_2_color_GREEN),
  #                 lm(range_RGB_color~range_2_color_BLUE))
  
  lm_RED<-lm(range_RGB_color~range_2_color_RED)
  lm_GREEN<-lm(range_RGB_color~range_2_color_GREEN)
  lm_BLUE<-lm(range_RGB_color~range_2_color_BLUE)
  
  #getting slope and intercept
  slope_all<-c(lm_RED$coefficients[2],
               lm_GREEN$coefficients[2],
               lm_BLUE$coefficients[2])
  
  intercept_all<-c(lm_RED$coefficients[1],
                   lm_GREEN$coefficients[1],
                   lm_BLUE$coefficients[1])
  
  #normalizing pca results to RGB for each Component (RED,GREEN, and BLUE)
  
  
  mRed<-(Main_3_Comps[,1]*slope_all[1])+intercept_all[1]
  
  mGreen<-(Main_3_Comps[,2]*slope_all[2])+intercept_all[2]
  
  mBlue<-(Main_3_Comps[,3]*slope_all[3])+intercept_all[3]
  
  new_Matriz<-cbind(mRed,mGreen,mBlue)
  
  #if there are negative values => put a zero in their positions
  
  new_Matriz[new_Matriz<0]<-0
  
  # Generating Hex Color Codes
  
  Hex_Color_Vector<-vector()
  
  row_Final<-length(new_Matriz[,1])
  
  columns_Final<-length(new_Matriz[1,])
  
  for(i in 1:row_Final){
    
    for(j in 1:columns_Final){
      
      red<-new_Matriz[i,1]
      green<-new_Matriz[i,2]
      blue<-new_Matriz[i,3]
      
      Hex_Color_Vector[i]<-rgb(red,green,blue,maxColorValue = 1)
      
    }
    
  }
  
  return(Hex_Color_Vector)
  
  
}

Creating_image<-function(nu_rows,nu_cols,len_for_color_matriz,Hex_color_vector){
  
  # loading library
  library(grid)
  
  # Creating image
  grid_matrix<-matrix(unlist(Hex_color_vector[1:len_for_color_matriz]), 
                      nrow=nu_rows, ncol = nu_cols, byrow = TRUE)
  
  
  imageFinal<-grid.raster(grid_matrix, interpolate=FALSE)
  
  return(imageFinal)
  
}

FUN <- function(x) {
  x <- as.integer(x)
  div <- seq_len(abs(x))
  factors <- div[x %% div == 0L]
  factors <- list(neg = -factors, pos = factors)
  return(factors)
}

############# MODIFIED FUNCTIONS ####################

Getting_first_three_components_MODI_prcomp<-function(matrix_for_PCA){
  
  # doing PCA with prcomp
  prcomp_All<-prcomp(matrix_for_PCA,
                     retx=TRUE)
  
  #obtaing the sum of varianza of the first three components
  
  Var<-prcomp_All$sdev^2 / sum(prcomp_All$sdev^2)
  
  Variance_all_3<-Var[1:3]
  
  Suma_total<-sum((Variance_all_3)*100)
  
  ## obteniendo coord of individuals of first three components 
  
  loadings_prcomp<-prcomp_All$x
  
  #   coordinates of individuals for the three colors
  
  Main_3_Comps<-loadings_prcomp[,c(1,2,3)]
  
  Variances_vs_dim_plot<-plot(prcomp_All,type="l")
  
  return(list("Var_C1_Red"=Variance_all_3[1],"Var_C2_Green"=Variance_all_3[2],
              "Var_C3_Blue"=Variance_all_3[3],"Total_First_3_Comps"=Suma_total,
              "Main_3_Comps"=Main_3_Comps,"pca"=prcomp_All,
              "Variances_vs_dim_plot"=Variances_vs_dim_plot))
  
}


Getting_first_three_components_MODI_PCA<-function(matrix_for_PCA){
  
  # looading library FactoMineR
  library(FactoMineR)
  
  #doing PCA with PCA function
  PCA_for_3_dimensions<-PCA(matrix_for_PCA,scale.unit=FALSE,
                            ncp = 3,graph=FALSE)
  
  # cumulative variance
  Suma_total_variance_3_comps<-PCA_for_3_dimensions$eig[,3][3]
  
  # coordinates of individuals for the three colors
  coordinates_of_individuals<-PCA_for_3_dimensions$ind$coord
  
  #return objects
  return(list("Var_C1_Red"=coordinates_of_individuals[,1],
              "Var_C2_Green"=coordinates_of_individuals[,2],
              "Var_C3_Blue"=coordinates_of_individuals[,3],
              "Total_First_3_Comps"=Suma_total_variance_3_comps,
              "Main_3_Comps"=coordinates_of_individuals,
              "PCA"=PCA_for_3_dimensions))
  
  
}

# making function

# default parameters : pedazos = 2000 porque se hara una imagen de 20 filas por 100 columnas
#                     sobrelapamiento = 2000 en relacion 1 a los pedazos 


Obtaining_window_to_slide_with_Overlapping<-function(genome_size)
{
  # pedazos es por defecto para 20 filas y 100 columnas
  pedazos<-2000
  
  window<-as.integer(genome_size/(pedazos*(1-0.5)))
  
  sobrelapamiento<-as.integer(window*0.5)
  
  return(list("window"=window,"sobrelapamiento"=sobrelapamiento))
  
}

Obtaining_window_to_slide_OVERLAPPING_cuadrada<-function(list_geno)
{
  # pedazos es por defecto para 25 filas y 25 columnas
  pedazos<-2025  
  
  window_list<-list()
  
  for(key in names(list_geno)){
    corte<-as.integer(list_geno[[key]]/pedazos)
    
    sobrelape<-as.integer(corte*0.5)
    
    window_list[[key]]<-c(corte,sobrelape)
  }
  
  return(window_list)
  
}

# Obtaining_window_to_slide_OVERLAPPING_cuadrada(archae_list)


Obtaining_window_to_slide_NO_OVERLAPPING_cuadrada<-function(list_genomes_size)
{
  pedazos<-2025
  
  window_list<-list()
  
  for(key in names(list_genomes_size)){
    window_list[[key]]<-as.integer(list_genomes_size[[key]]/pedazos)
  }
  return(window_list)
}

############################################ FUNCTIONS #######################################################


#1 data genome should be transformed to data.frame object

Creating_density_plots_for_each_nucleotide<-function(data_frame_from_genome,OligoNames)
{
  # para A 1-64 ; C 65-128 ; G 129-192 ; T 193-256
  
  iniciales<-vector()
  finales<-vector()
  
  for(i in 1:16)
  {
    final=16*i
    
    inicial = final-15
    
    if(final%%i==0)
    {
      
      # salta de a 4 posiciones en inicial y finales para cada nucleotido
      iniciales[i]<-inicial
      finales[i]<-final
      #         writeLines(paste("inicial =", inicial, " final=", final))
    }
  }
  
  # making data for density plot for all nucleotides
  
  ###### FOR A  ######
  
  data_AA_all_f<-data_frame_from_genome[iniciales[1]:finales[1]]
  AA_all_f<-OligoNames[[1]][iniciales[1]:finales[1]]
  
  data_AC_all_f<-data_frame_from_genome[iniciales[2]:finales[2]]
  AC_all_f<-OligoNames[[1]][iniciales[2]:finales[2]]
  
  data_AG_all_f<-data_frame_from_genome[iniciales[3]:finales[3]]
  AG_all_f<-OligoNames[[1]][iniciales[3]:finales[3]]
  
  data_AT_all_f<-data_frame_from_genome[iniciales[4]:finales[4]]
  AT_all_f<-OligoNames[[1]][iniciales[4]:finales[4]]
  
  data_for_density_plot_A<-cbind(data_AA_all_f,data_AC_all_f,data_AG_all_f,
                                 data_AT_all_f)
  names_for_density_plot_A<-c(AA_all_f,AC_all_f,AG_all_f,AT_all_f)
  
  data_for_A<-list(data_for_density_plot_A,names_for_density_plot_A)
  
  ###### FOR C  ######
  
  data_CA_all<-data_frame_from_genome[iniciales[5]:finales[5]]
  CA_all<-OligoNames[[1]][iniciales[5]:finales[5]]
  
  data_CC_all<-data_frame_from_genome[iniciales[6]:finales[6]]
  CC_all<-OligoNames[[1]][iniciales[6]:finales[6]]
  
  data_CG_all<-data_frame_from_genome[iniciales[7]:finales[7]]
  CG_all<-OligoNames[[1]][iniciales[7]:finales[7]]
  
  data_CT_all<-data_frame_from_genome[iniciales[8]:finales[8]]
  CT_all<-OligoNames[[1]][iniciales[8]:finales[8]]
  
  data_for_density_plot_C<-cbind(data_CA_all,data_CC_all,data_CG_all,data_CT_all)
  names_for_density_plot_C<-c(CA_all,CC_all,CG_all,CT_all)
  
  data_for_C<-list(data_for_density_plot_C,names_for_density_plot_C)
  
  ###### FOR G  ######
  
  data_GA_all<-data_frame_from_genome[iniciales[9]:finales[9]]
  GA_all<-OligoNames[[1]][iniciales[9]:finales[9]]
  
  data_GC_all<-data_frame_from_genome[iniciales[10]:finales[10]]
  GC_all<-OligoNames[[1]][iniciales[10]:finales[10]]
  
  data_GG_all<-data_frame_from_genome[iniciales[11]:finales[11]]
  GG_all<-OligoNames[[1]][iniciales[11]:finales[11]]
  
  data_GT_all<-data_frame_from_genome[iniciales[12]:finales[12]]
  GT_all<-OligoNames[[1]][iniciales[12]:finales[12]]
  
  data_for_density_plot_G<-cbind(data_GA_all,data_GC_all,data_GG_all,data_GT_all)
  names_for_density_plot_G<-c(GA_all,GC_all,GG_all,GT_all)
  
  data_for_G<-list(data_for_density_plot_G,names_for_density_plot_G)
  
  ###### FOR T  ######
  
  data_TA_all<-data_frame_from_genome[iniciales[13]:finales[13]]
  TA_all<-OligoNames[[1]][iniciales[13]:finales[13]]
  
  data_TC_all<-data_frame_from_genome[iniciales[14]:finales[14]]
  TC_all<-OligoNames[[1]][iniciales[14]:finales[14]]
  
  data_TG_all<-data_frame_from_genome[iniciales[15]:finales[15]]
  TG_all<-OligoNames[[1]][iniciales[15]:finales[15]]
  
  data_TT_all<-data_frame_from_genome[iniciales[16]:finales[16]]
  TT_all<-OligoNames[[1]][iniciales[16]:finales[16]]
  
  data_for_density_plot_T<-cbind(data_TA_all,data_TC_all,data_TG_all,data_TT_all)
  names_for_density_plot_T<-c(TA_all,TC_all,TG_all,TT_all)
  
  data_for_T<-list(data_for_density_plot_T,names_for_density_plot_T)
  
  ###### returning all elements ######
  
  inicialesForDensityPlot<-iniciales[1:4]
  finalesForDensityPlot<-finales[1:4]
  
  return(list("data_for_A"=data_for_A,
              "data_for_C"=data_for_C,
              "data_for_G"=data_for_G,
              "data_for_T"=data_for_T
  ))
  
}

######### FUNCION Building_density_plot

Building_density_plot<-function(nameN,data_from_nucleotide,cexT)
{
  main_names<-c(paste(nameN,"A*",sep=""),paste(nameN,"C*",sep=""),paste(nameN,"G*",sep=""),paste(nameN,"T*",sep=""))
  
  ancho=NULL
  largo=NULL
  
  ## iniciales y finales
  
  iniciales<-c(1,17,33,49)
  finales<-c(16,32,48,64)
  
  par(mfrow=c(1,4))
  
  for(i in 1:4)
  {
    
    for(j in 1:length(data_from_nucleotide[[1]][iniciales[i]:finales[i]]))
    {
      ancho = c(ancho, 
                density(data_from_nucleotide[[1]][iniciales[i]:finales[i]][[j]])$x)
      largo = c(largo, 
                density(data_from_nucleotide[[1]][iniciales[i]:finales[i]][[j]])$y)
    }
    
    anchoR<-range(ancho)
    largoR<-range(largo)
    
    # all(data_for_Density_plots$data_for_A[[1]][1]==data_AA_all[[1]])
    
    plot(density(data_from_nucleotide[[1]][iniciales[i]:finales[i]][[j]]),
         xlim = anchoR, ylim = largoR, main = main_names[i])
    
    c16<-rainbow(22)
    legendCol<-vector()
    
    for(j in 1:length(data_from_nucleotide[[1]][iniciales[i]:finales[i]]) ){
      lines(density(data_from_nucleotide[[1]][iniciales[i]:finales[i]][[j]]), 
            xlim = anchoR, ylim = largoR, col = c16[j])
      legendCol[j]<-c16[j]
    }
    
    # Add legend to top right, outside plot region
    legend("topright", inset=c(-0.001,0), 
           legend=c(data_from_nucleotide[[2]][iniciales[i]:finales[i]]),fill= legendCol,cex = cexT)
    
  }
  
}

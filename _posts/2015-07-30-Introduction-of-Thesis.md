---
layout: post
title: Introduction of Thesis
description: "Introduction"
category: Thesis
tags: [Introduction]
---
{% include JB/setup %}

# Intro

------

The objective of my thesis is to create an image from a genome but this image should have a biological insight. Because humans are more capable of finding patterns visually than by reading a plain text, creating an image from a genome is a great way of understanding the data. Also, there is not a standardized methodology on the literature related with my project, so my thesis director and I made a draft of a methodology that we want to try as specified in this file [preliminary project in Spanish](/Additional_material/Preliminary_Project_Camila_Martinez.docx)

------

## Objective of my thesis

**_Main objective:_** Create a function that takes as input **compositional biases of genomes** and throws as output an image with biological sense.

------

## What are compositional biases of genomes

A genome is a sequence of characters with an **alfabeth = {A,T,G,C}** so it is possible to count the frequency of each element or a combination of these elements. Then, it is important the length of the sequence that is going to be counted, which is denominated as **k-mer**. 

Also according to publications such as **_Compositional biases of bacterial genomes and evolutionary implications by S  Karlin, J Mrázek and A M Campbel in 1997_** and **_Evolutionary Implications of Microbial Genome Tetranucleotide Frequency Biases by David T. Pride, Richard J. Meinersmann, Trudy M. Wassenaar and Martin J. Blase in 2003_**, the different k-mer frequency biases have evolutionary consequences that can describe differences between organisms without resorting to methods based on homology comparison. 

------

## Components of the main function

#### Note: I have decided to divide my main function in 6 smaller function that make different tasks:

1. Function1: Processing_fragments.
     
      * Input <= **_PATH_** where the fragmented genome of a organism is.
      
      * Output => **_List_** that contains all the fragments of a genome.

2. Function2: GettingFrequency_account.
     
      * Input <= List of all the fragments, and the **_k-mer_** that is going to be counted.
      
      * Output => Matrix where the **_rows_** are the number of the fragments and the **_columns_** are the possible combinations of the chosen k-mer. For example, if the k-mer = 4, then there are 256 possible combinations (4^4 = 256).

3. Function3: GettingFirstThreeComponents_MODI.
     
      * Input <= The matrix that resulted from the function before.
      
      * Output => R list that has different objects:

            3A. The variance of the three components resulted from the PCA, where component1 = Red, component2 = Green, and component3 = Blue.
            3B. The total sum of the variances of the first 3 components.
            3C. A matrix, where the rows are the names of the fragments and the columns are the first three components for each fragment.

4. Function4: GeneratingHexVector_MODI.
     
      * Input <= A matrix that has the first three components of the PCA for each fragment of the initial genome
      
      * Output => A vector of **_length equals to the number of fragments_**, where each element has the hex color code of each fragment

5. Function5 and Function6: ObtainingDimensionsForImage(HexColorVector) and Creating_Image(nRows,nCols,lenColorMatrix,HexColorVector)
     
      * These two functions work to create the image.














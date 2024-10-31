# Deciphering progressive lesion areas in breast cancer spatial transcriptomics via TGR-NMF
# Overview
![flowchat](https://github.com/xiangshanxs/TGR-NMF/blob/main/TGR-NMF.jpg)
TGR-NMF is an interpretable model for identifying interesting spatial domains that reveal breast cancer progression by using spatial transcriptomic data. TGR-NMF accounts for the impact of dropout events in computation and introduces three graph regularization terms by integrating gene expression count data and positional coordinate information. The workflow for deciphering spatial domains in breast cancer spatial transcriptomics is shown in fig. Firstly, the unitization strategy is adopted to process the gene expression count matrix and position coordinate matrix. Secondly, TGR-NMF is constructed by integrating three distinct graph regularization terms to obtain a low-dimensional representation of the gene expression count data. The construction of these regularization terms involves heat kernel weighting for the gene expression count matrix, cosine similarity and heat kernel weighting for the position matrix. Next, the K-means algorithm is employed to identify the spatial domains based on the obtained low-dimensional representation. Finally, progressive lesion regions of breast cancer are revealed through various downstream analyses.
# Datasets
All datasets used in our paper can be found in:
The first HBRC spatial transcriptomic dataset is collected from the 10x Genomics website at https://www.10xgenomics.com/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0. The gold standard of the first HBRC dataset is accessible at https://github.com/JinmiaoChenLab/SEDR_analyses/blob/master/data/BRCA1/metadata.tsv. The second HBRC spatial transcriptomic dataset is downloaded from https://zenodo.org/records/4739739.
# Usage
We provided some demos to demonstrate usage of TGR-NMF.

    clc
    clear;
    close all;  
    rng(0)
    idx = randperm(3798);
    load('C:/Users/Administrator/Desktop/BRCA_exp.mat') 
    load('C:/Users/Administrator/Desktop/BRCA_pos2.mat') 
    X1 = Unitization(in_X);   
    X2 = Unitization(pos);
    nnClass = 20;     
    GSCA = [];
    NMI = [];
%% MGNMF
%============================= Construct W1 ===============================
        options = [];
        options.NeighborMode = 'KNN';  
        options.k = 30;    
        options.WeightMode = 'HeatKernel'; 
        options.t = 1;     
    W1 = ConstructW1(X1,options); 
%============================ Construct W2 ===============================
        options = [];
        options.NeighborMode = 'KNN';   
        options.k = 9;     
        options.WeightMode = 'Cosine'; 
    W2 = ConstructW2(X2,options); 
%============================= Construct W3 ===============================
        options = [];
        options.NeighborMode = 'KNN';  
        options.k = 9;     
        options.WeightMode = 'HeatKernel'; 
        options.t = 1;    
    W3 = ConstructW3(X2,options);  

      [mFea,nSmp]=size(X1');
      k = 1190;
      U = abs(rand(mFea,k));    
      V = abs(rand(k,nSmp));    
        options = [];
        options.alpha = 40;
        options.mu1 = 0.8;    
        options.mu2 = 0.01;    
        options.mu3 = 0.19;    
        options.error = 0.005;   
        options.maxIter = 500;    
        options.nRepeat = 1;    
        options.minIter = 20;  
        options.meanFitRatio = 0.5;   
        [U2, V2] = GNMF_Multi(X1', k, W1, W2, W3, options, U, V');
        [project_labs, center] = Litekmeans(V2, nnClass, 'Start', idx(1:nnClass)); 
%============================= evaluate ==================================
        ARI = Cal_ARI(true_labs, project_labs);
        GSCA = GSCA_ClusteringMeasure(true_labs, project_labs);
        

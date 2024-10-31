    clc
    clear;
    close all;  
    rng(0)
    idx = randperm(3798);
    load('C:/Users/Administrator/Desktop/in_X.mat') 
    load('C:/Users/Administrator/Desktop/pos.mat') 
    load('C:/Users/Administrator/Desktop/true_labs.mat')
    X1 = Unitization(in_X);   
    X2 = Unitization(pos);
    nnClass = 20;     
    GSCA = [];
    NMI = [];
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
        

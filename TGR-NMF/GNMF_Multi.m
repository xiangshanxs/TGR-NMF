function [U_final, V_final, nIter_final, objhistory_final] = GNMF_Multi(X, k, W1, W2, W3, options, U, V)  
differror = options.error;  
maxIter = options.maxIter;
nRepeat = options.nRepeat;  
minIter = options.minIter - 1;   
meanFitRatio = options.meanFitRatio; 
alpha = options.alpha;  
mu1 = options.mu1;
mu2 = options.mu2;
mu3 = options.mu3;
Norm = 2;  
NormV = 0;  
[mFea,nSmp]=size(X);    

    %======================================================
    W1 = alpha*mu1*W1;
    DCol = full(sum(W1,2));  
    D1 = spdiags(DCol,0,nSmp,nSmp);  
    L1 = D1 - W1;   
    %======================================================
    W2 = alpha*mu2*W2;
    DCol = full(sum(W1,2));  
    D2 = spdiags(DCol,0,nSmp,nSmp);  
    L2 = D2 - W2;   
    %======================================================
    W3 = alpha*mu3*W3;
    DCol = full(sum(W3,2));  
    D3 = spdiags(DCol,0,nSmp,nSmp);  
    L3 = D3 - W3;   
    nRepeat = 1;    
[U,V] = UnitizationUV(U, V, NormV, Norm);   
    selectInit = 0;  
    minIter = 0;
    tryNo = 0;    
    nIter = 0;   
while tryNo < nRepeat  
    tryNo = tryNo+1;
    maxErr = 1;
    while(maxErr > differror)  
      % ===================== update U ========================
        XV = X*V;   
        VV = V'*V; 
        UVV = U*VV; 
        
        U = U.*(XV./max(UVV,1e-10)); % 3mk
      % ===================== update V ========================
        XU = X'*U;  
        UU = U'*U;  
        VUU = V*UU; 
        
        if alpha > 0
            W1V = W1*V;
            W2V = W2*V;
            W3V = W3*V;
            D1V = D1*V;
            D2V = D2*V;
            D3V = D3*V;

            XU = XU + W1V + W2V + W3V;       
            VUU = VUU + D1V + D2V + D3V;     
        end
        V = V.*(XU./max(VUU,1e-10)); 
        
      % ===================== update  ========================       
        nIter = nIter + 1;
        if nIter > minIter   
            maxErr = 1;  
            if nIter >= maxIter  
               maxErr = 0;  
               objhistory = 0;  
            end
        end
    end
    
    if tryNo == 1   
        U_final = U;  
        V_final = V;
        nIter_final = nIter;  
        objhistory_final = objhistory;
    else 
       if objhistory(end) < objhistory_final(end)  
           U_final = U;
           V_final = V;
           nIter_final = nIter;
           objhistory_final = objhistory;
       end
    end
end

[U_final,V_final] = UnitizationUV(U_final, V_final, NormV, Norm);


function [U, V] = UnitizationUV(U, V, NormV, Norm)
    K = size(U,2);
    if Norm == 2
            norms = max(1e-15,sqrt(sum(U.^2,1)))';   
            U = U*spdiags(norms.^-1,0,K,K); 
            V = V*spdiags(norms,0,K,K);  
     end

        

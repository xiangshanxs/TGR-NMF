function W = ConstructW2(fea,options)  
%=================================================
   if ~isfield(options,'bSelfConnected')
        options.bSelfConnected = 0;   
   end
   if ~isfield(options,'bUnitizated')
    options.bUnitizated = 0;
   end   
%=================================================
      bCosine = 1;
      nSmp = size(fea,1);  
   if bCosine && ~options.bUnitizated   
      Unitfea = Unitization(fea);
   end
%======================= ==========================
        G = zeros(nSmp*(options.k+1),3);  
        smpIdx = 1:nSmp;   
        dist = EuDist2(fea,fea,0);     
        [dump, idx] = sort(dist,2); 
        idx = idx(:,1:options.k+1);  
        dump = dump(:,1:options.k+1); 
%==================================================
    if bCosine
       dist = Unitfea*Unitfea'; 
       dist = full(dist);
       linidx = (1:nSmp)'; % 3798×1
       dump = dist(sub2ind(size(dist),linidx(:,ones(1,size(idx,2))),idx));  
    end

      G(1:nSmp*(options.k+1),1) = repmat(smpIdx',[options.k+1,1]);
      G(1:nSmp*(options.k+1),2) = idx(:);
      G(1:nSmp*(options.k+1),3) = dump(:);
   W = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
    if ~options.bSelfConnected
       W = W - diag(diag(W));
    end
        W = max(W,W');
return;    

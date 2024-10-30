function W = ConstructW3(fea,options)  
%=================================================
   if ~isfield(options,'bSelfConnected')
        options.bSelfConnected = 0;   
   end
%=================================================
        nSmp = size(fea,1);  
        G = zeros(nSmp*(options.k+1),3);  
        smpIdx = 1:nSmp;   
        dist = EuDist2(fea,fea,0);     
        dump = zeros(nSmp,options.k+1);  
        idx = dump;    
   for j = 1:options.k+1  
         [dump(:,j),idx(:,j)] = min(dist,[],2);  
         temp = (idx(:,j)-1)*nSmp+(1:nSmp)';  
         dist(temp) = 1e100;  
   end    
          dump = exp(-dump/(2*options.t^2));
          G(1:nSmp*(options.k+1),1) = repmat(smpIdx',[options.k+1,1]);
          G(1:nSmp*(options.k+1),2) = idx(:);
          G(1:nSmp*(options.k+1),3) = dump(:);
          W = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);  
    if ~options.bSelfConnected
        W = W - diag(diag(W));   
    end
          W = max(W,W');
          return;

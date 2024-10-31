function fea = Unitization(fea,row,norm)

if ~exist('row','var')      
    row = 1;
end

if ~exist('norm','var')     
    norm = 2;
end
if norm < 1
    error('It is not a norm when p small than 1!');  
end 
    nSmp = size(fea,1);   
    feaNorm = max(1e-14,full(sum(abs(fea).^norm,2)));       
    fea = spdiags(feaNorm.^-(1./norm),0,nSmp,nSmp)*fea;          
return;
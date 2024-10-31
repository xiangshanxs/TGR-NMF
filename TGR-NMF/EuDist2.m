function D = EuDist2(fea_a,fea_b,bSqrt)

    aa = sum(fea_a.*fea_a,2); %a×1
    bb = sum(fea_b.*fea_b,2); %b×1
    ab = fea_a*fea_b'; %a×b

    if issparse(aa)
        aa = full(aa);
        bb = full(bb);
    end

    D = bsxfun(@plus,aa,bb') - 2*ab;
    D(D<0) = 0;

function [Z,A]=genZA(G,maf)
    [m,n]=size(G);
    Z=zscore_sv(G,1,-1,maf); %zscore by maf excluding missing.. then replace missing with zero
    a=Z*Z'/n;
    %adjust diagonal by inbreeding effect, for non-gen appplications A=a
    v=(2*maf-1)./sqrt(2*maf.*(1-maf)); 
    Q=Z*v;
    Q=Q/n;
    Q=eye(m).*repmat(Q,[1,m]);
    A=a+Q;
end
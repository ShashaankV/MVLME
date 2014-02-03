function X=genX(Ntr,C)
    [m,n]=size(C);
    X=zeros(m*Ntr,n*Ntr);
    for i=1:Ntr
        indm=(i-1)*m+1:(i-1)*m+m;
        indn=(i-1)*n+1:(i-1)*n+n;
        X(indm,indn)=C;
    end
end


    

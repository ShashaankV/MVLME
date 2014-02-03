
function [Finv,Finvmap]=FisherInv(A,X,covg,cove)
    ntraits = size(covg,1);
    m=size(A,1);
    I=eye(m,'single');
    
    A=single(A);

    Vinv=zeros(m*ntraits,'single');
    Vtinv=zeros(m*ntraits,'single');
    BigQ=zeros(m*ntraits,'single');
    

    if ntraits > 1
        
    [S,cgt]=eig((covg/cove)');

    Q=sqrt(S'*cove*S)\S';
    
    
    for i=1:ntraits
        i_ind=(i-1)*m+1:i*m;
        Vtinv(i_ind,i_ind)=(cgt(i,i)*A+I)\I;
    end


    for k=1:ntraits
        k_ind=(k-1)*m+1:k*m;
        for i=1:ntraits
        i_ind=(i-1)*m+1:i*m;
           for j=1:ntraits
            j_ind=(j-1)*m+1:j*m;
                Vinv(i_ind,j_ind)=Vinv(i_ind,j_ind)+Q(k,i)*Vtinv(k_ind,k_ind)*Q(k,j);
            end
        end
    end
    
    Vinv=real(Vinv);
    
    clear Vtinv
    else
        V=covg*A+cove*I;
        Vinv=I/V;
    end
      
    VX=Vinv*X;
    XVX=X'*VX;
    P=Vinv-VX/XVX*VX';
    
    P=single(P);
    
    clear Vinv X VX XVX
    

        for i=1:ntraits
        i_ind=(i-1)*m+1:i*m;
        for j=i:ntraits
            j_ind=(j-1)*m+1:j*m;
            for k=1:ntraits
                k_ind=(k-1)*m+1:k*m;
                for l=k:ntraits
                    l_ind=(l-1)*m+1:l*m;
                    VAij=zeros(size(P),'single');
                    VEij=zeros(size(P),'single');
                    VAkl=zeros(size(P),'single');
                    VEkl=zeros(size(P),'single');
                    VAij(i_ind,j_ind)=A;
                    VEij(i_ind,j_ind)=I;
                    if i~=j
                        VAij(j_ind,i_ind)=A;
                        VEij(j_ind,i_ind)=I;
                    end
                    VAkl(k_ind,l_ind)=A;
                    VEkl(k_ind,l_ind)=I;
                    if k~=l
                        VAkl(l_ind,k_ind)=A;
                        VEkl(l_ind,k_ind)=I;
                    end
                                        
                    a=ntraits*(i-1)-(i>1)*(i-1)*(i-2)/2+j-i+1;
                    b=ntraits*(k-1)-(k>1)*(k-1)*(k-2)/2+l-k+1;
                    
                    PVAij=P*VAij;
                    PVAkl=P*VAkl;
                    PVEij=P*VEij;
                    PVEkl=P*VEkl;
                    
                     
                    F(2*a-1,2*b-1)=trace(PVAij*PVAkl);
                    Finvmap(2*a-1,2*b-1)={['g' num2str(i) 'g' num2str(j) '_g' num2str(k) 'g' num2str(l)]};
                    
                    F(2*a-1,2*b)=trace(PVAij*PVEkl);
                    Finvmap(2*a-1,2*b)={['g' num2str(i) 'g' num2str(j) '_e' num2str(k) 'e' num2str(l)]};
                    
                    F(2*a,2*b-1)=trace(PVEij*PVAkl);
                    Finvmap(2*a,2*b-1)={['e' num2str(i) 'e' num2str(j) '_g' num2str(k) 'g' num2str(l)]};
                    
                    F(2*a,2*b)=trace(PVEij*PVEkl);
                    Finvmap(2*a,2*b)={['e' num2str(i) 'e' num2str(j) '_e' num2str(k) 'e' num2str(l)]};
                    
                end
            end
        end
    end
    

    F=F/2;
    Finv=inv(F);
    

end


function [dLL,covg,cove]=REMLEM(A,y,X,covg,cove,alpha,varargin)
    if nargin<6
        alpha=1;
    end

    ntraits = size(covg,1);
    m=size(A,1);
    I=eye(m,'single');
    


    A=single(A);

    Vinv=zeros(m*ntraits,'single');
    Vtinv=zeros(m*ntraits,'single');
    BigQ=zeros(m*ntraits,'single');
    

    [S,cgt]=eig((covg/cove)');

    Q=sqrt(S'*cove*S)\S';
    
    
    for i=1:ntraits
        i_ind=(i-1)*m+1:i*m;
        Vtinv(i_ind,i_ind)=(cgt(i,i)*A+I)\I;
    end
    clear I


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
      
    
    
    VX=Vinv*X;
    XVX=X'*VX;
    P=Vinv-VX/XVX*VX';
    
    clear Vinv;


    %update parameters
       dLL=0;
    for i=1:ntraits
      i_ind=(i-1)*m+1:i*m;
        for j=i:ntraits
            j_ind=(j-1)*m+1:j*m;
            Bg=0;
            Be=0;
            for k=1:ntraits
                k_ind=(k-1)*m+1:k*m;
                for l=1:ntraits    
                    l_ind=(l-1)*m+1:l*m;
                    ui=y(:,k)'*P(k_ind,i_ind);
                    vj=P(j_ind,l_ind)*y(:,l);
                    Bg=Bg+ui*A*vj;
                    Be=Be+ui*vj;
                end
            end
            dLg=-trace(P(i_ind,j_ind)*A)+Bg;
            dLe=-trace(P(i_ind,j_ind))+Be;
            
            covg(i,j)=covg(i,j)*(1+alpha*covg(i,j)/m*dLg);
            cove(i,j)=cove(i,j)*(1+alpha*cove(i,j)/m*dLe);
                                
            
            LLg=(i~=j+1)*covg(i,j)^2/m*dLg^2;       
            LLe=(i~=j+1)*cove(i,j)^2/m*dLe^2;
            
            dLL=dLL+LLg+LLe;
        end
    end
    
    covg=real(triu(covg)+triu(covg,1)');
    cove=real(triu(cove)+triu(cove,1)');
    
    clear P ui vj

    
end


function [dLL,covg,cove]=REMLEM_trans(A,y,X,covg,cove)


ntraits=size(covg,1);
nfe=size(X,2)/ntraits; %number of fixed effects
m=size(A,1);

I=eye(m);

% canonical transformation
[S,cgt]=eig((covg/cove)');

Q=sqrt(S'*cove*S)\S';

yv=y*Q';
dLL=0;
for i=1:ntraits
    %grab fixed effects block for this trait
    i_ind=(i-1)*m+1:i*m; 
    j_ind=(i-1)*nfe+1:i*nfe;
    x=X(i_ind,j_ind);
    
    V=cgt(i,i)*A+I;
    Vinv=V\I;

    Vx=Vinv*x;
    xVx=x'*Vx;
    P=Vinv-Vx/xVx*Vx';
    
    clear Vinv;
    
    PVg=P*A;
    
    % Bg and Be are zero at maximum likelihood
    
    Bg=-trace(PVg)+yv(:,i)'*PVg*P*yv(:,i);
    Be=-trace(P)+yv(:,i)'*P*P*yv(:,i);
    
    
    cgt(i,i)=cgt(i,i)*(1 + cgt(i,i)*Bg/m);
    cet(i,i)=(1 + Be/m);
    
    LLg=abs(cgt(i,i)^2/m*Bg^2);       
    LLe=abs(cet(i,i)^2/m*Be^2);
            
    dLL=dLL+LLg+LLe;
            
    
            
end
    

covg=real(Q\cgt/Q');
cove=real(Q\cet/Q');


end

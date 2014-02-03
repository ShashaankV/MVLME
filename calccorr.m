function [r,rse,rpval]=calccorr(tr1,tr2,eclass,Sigma,Finv,Finvmap)
    %grab relevant parameter estimates from Sigma
    v1=Sigma(tr1,tr1);
    v2=Sigma(tr2,tr2);
    cov12=Sigma(tr1,tr2);
    %grab relevant variance/covariance around the estimated means from Finv
    %to do this create the element label, find it in Fmap, and grab it from
    %Finv
    %1. variance of v1
    lab=[eclass num2str(tr1) eclass num2str(tr1) '_' eclass num2str(tr1) eclass num2str(tr1)];
    ind=find(strcmp(Finvmap,lab)==1);
    var_v1=Finv(ind);
    
    lab=[eclass num2str(tr2) eclass num2str(tr2) '_' eclass num2str(tr2) eclass num2str(tr2)];
    ind=find(strcmp(Finvmap,lab)==1);
    var_v2=Finv(ind);
    
    lab=[eclass num2str(tr1) eclass num2str(tr2) '_' eclass num2str(tr1) eclass num2str(tr2)];
    ind=find(strcmp(Finvmap,lab)==1);
    var_cov12=Finv(ind);
    
    lab=[eclass num2str(tr1) eclass num2str(tr1) '_' eclass num2str(tr2) eclass num2str(tr2)];
    ind=find(strcmp(Finvmap,lab)==1);
    cov_v1v2=Finv(ind);
    
    lab=[eclass num2str(tr1) eclass num2str(tr1) '_' eclass num2str(tr1) eclass num2str(tr2)];
    ind=find(strcmp(Finvmap,lab)==1);
    cov_v1cov12=Finv(ind);
    
    lab=[eclass num2str(tr2) eclass num2str(tr2) '_' eclass num2str(tr1) eclass num2str(tr2)];
    ind=find(strcmp(Finvmap,lab)==1);
    cov_v2cov12=Finv(ind);
    
    r=cov12/sqrt(v1*v2);
    
    rse=sqrt(r^2*(var_v1/(4*v1^2)+var_v2/(4*v2^2)+var_cov12/(cov12^2)+2*cov_v1v2/(4*v1*v2)-2*cov_v1cov12/(2*v1*cov12)-2*cov_v2cov12/(2*v2*cov12)));
    
    if sign(r)==-1
        rtmp=-r;
    else
        rtmp=r;
    end
    rpval=2*normcdf(-atanh(rtmp)/atanh(rse),0,1);
end

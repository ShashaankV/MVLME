function [h2,h2se]=calch2(tr,Sigma_g,Sigma_e,Finv,Finvmap)
    
    lab=['g' num2str(tr) 'g' num2str(tr) '_' 'g' num2str(tr) 'g' num2str(tr)];
    ind=find(strcmp(Finvmap,lab)==1);
    var_vg=Finv(ind);
    
    lab=['e' num2str(tr) 'e' num2str(tr) '_' 'e' num2str(tr) 'e' num2str(tr)];
    ind=find(strcmp(Finvmap,lab)==1);
    var_ve=Finv(ind);
    
    lab=['g' num2str(tr) 'g' num2str(tr) '_' 'e' num2str(tr) 'e' num2str(tr)];
    ind=find(strcmp(Finvmap,lab)==1);
    cov_vgve=Finv(ind);
    
    cov_vgvp=var_vg+cov_vgve; % the cov(vg,vp) = var(vg) + cov(vg,ve)
    
    vg=Sigma_g(tr,tr);
    ve=Sigma_e(tr,tr);
    vp=vg+ve;
    h2=vg/vp;
    
    var_vp=var_vg+var_ve+2*cov_vgve;
    
    var_h2=h2^2*(var_vg/vg^2+var_vp/vp^2-2*cov_vgvp/(vg*vp));
    
    h2se=sqrt(var_h2);
end
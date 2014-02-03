clc
clear
figure(1)
clf

load('simdat1')
Ntr=size(Y,2);
m=500; %to save on computation, only analyze the first m subjects (~500 is the min for simdat1)
G=G(1:m,:);
Y=zscore_sv(Y(1:m,:));
[Z,A]=genZA(G,maf);
% 
% 
%generate X 
COVARS=ones(m,1); %in this case we assume no covariates
X=genX(Ntr,COVARS);
% 
% 

%random initial covariance matrices
offdiag=(1-eye(Ntr)).*(rand(Ntr)-.5); 
offdiag=triu(offdiag)+triu(offdiag,1)';
Sigma_g=repmat(rand([Ntr,1]),[1 Ntr]).*eye(Ntr);
Sigma_e=(1-Sigma_g).*eye(Ntr);
Sigma_g=Sigma_g +offdiag; 
offdiag=(1-eye(Ntr)).*(rand(Ntr)-.5); 
offdiag=triu(offdiag)+triu(offdiag,1)'; 
Sigma_e=Sigma_e +offdiag;
 

COVARS=ones(m,1); %in this case we assume no covariates
X1=genX(1,COVARS);

tol=1*10^-4; 
disp('calculating univariate variances')
for tr=1:Ntr
    vg=rand;
    ve=rand;
    dl=1;
    while dl>tol
        [dl,vg,ve]=REMLEM(A,Y(:,tr),X1,vg,ve); 
    end
    Sigma_g(tr,tr)=vg;
    Sigma_e(tr,tr)=ve;
end





disp('initializing covariances with REMLEMtrans') 
i=1;
sigmagrec(:,:,i)=Sigma_g;
sigmaerec(:,:,i)=Sigma_e;

dl=1;
tol=1*10^-4;
while dl>tol
    i=i+1;
    [dl,Sigma_g,Sigma_e]=REMLEM_trans(A,Y,X,Sigma_g,Sigma_e);
    sigmagrec(:,:,i)=Sigma_g;
    sigmaerec(:,:,i)=Sigma_e;
end
itrans=i; %for plotting later

%use variances from univariate REMLEM, as REMLEMtrans may be off
for tr=1:Ntr
    Sigma_g(tr,tr)=squeeze(sigmagrec(tr,tr,1));
    Sigma_e(tr,tr)=squeeze(sigmaerec(tr,tr,1));
end

%main REML EM 

tol=1*10^-4; 
dl=1;
% 
'Starting REMLEM'
while dl>tol
    i=i+1;
    ['REMLEM dl: ' num2str(dl)] %this  can be commented out to save comp. time
    [dl,Sigma_g,Sigma_e]=REMLEM(A,Y,X,Sigma_g,Sigma_e); 
    sigmagrec(:,:,i)=Sigma_g;
    sigmaerec(:,:,i)=Sigma_e;
end

%plot search paths for the first two traits
%black vertical line is the REMLEMtrans, REMLEM boundary 
%h2 tr 1
vg1=squeeze(sigmagrec(1,1,:));
ve1=squeeze(sigmaerec(1,1,:));
h21_hat=vg1./(vg1+ve1);
%h2 tr 2
vg2=squeeze(sigmagrec(2,2,:));
ve2=squeeze(sigmaerec(2,2,:));
h22_hat=vg2./(vg2+ve2);
%rg tr 1 and 2
covg=squeeze(sigmagrec(1,2,:));
rg_hat=covg./sqrt(vg1.*vg2);
%rg tr 1 and 2
cove=squeeze(sigmaerec(1,2,:));
re_hat=cove./sqrt(ve1.*ve2);

figure(1)
clf
subplot(1,2,1)
plot(h21_hat,'linewidth',2)
hold on
plot(h22_hat,'r','linewidth',2)
axis([1 length(vg1) 0 1])
xlabel('iteration')
ylabel('h^2')
line([itrans,itrans],[0,1],'color','k')
h=legend({'$\hat{h}^2_{tr1}$','$\hat{h}^2_{tr2}$'},'fontsize',14,'location','Northwest');
set(h,'Interpreter','latex')

subplot(1,2,2)
plot(rg_hat,'g','linewidth',2)
hold on
plot(re_hat,'k','linewidth',2)
axis([1 length(vg1) -1 1])
xlabel('iteration')
ylabel('correlation')
line([itrans,itrans],[-1,1],'color','k')
h=legend({'$\hat{r}_g$','$\hat{r}_e$'},'fontsize',14,'location','Northwest');
set(h,'Interpreter','latex')
drawnow

% 
% 
%calc error
[Finv,Finvmap]=FisherInv(A,X,Sigma_g,Sigma_e);

%est h2  
for i=1:Ntr
    [H2hat(i),H2se(i)]=calch2(i,Sigma_g,Sigma_e,Finv,Finvmap);
end

%est C 
Cghat=eye(Ntr);
Cgse=zeros(Ntr);
Cgpval=zeros(Ntr);

Cehat=eye(Ntr);
Cese=zeros(Ntr);
Cepval=zeros(Ntr);

for i=1:Ntr
    for j=i+1:Ntr
        [Cghat(i,j),Cgse(i,j),Cgpval(i,j)]=calccorr(i,j,'g',Sigma_g,Finv,Finvmap);
        [Cehat(i,j),Cese(i,j),Cepval(i,j)]=calccorr(i,j,'e',Sigma_e,Finv,Finvmap);
    end
end

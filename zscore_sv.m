%zscore with missing values, default over rows
%1. calculate mean and std only over non-missing
%2. normalize only non-missing
%3. replace missing with zero
function A=zscore_sv(A,dim,miss_val,maf,varargin)
    Narg=nargin;
    if Narg>1
        if dim==2 %if over columns, then transpose input and transpose back before return
            A=A';
        end
    end
    [m,n]=size(A);
    for i=1:n
        if Narg>2 %if missing elements; locations given by miss_val
            ind=A(:,i)~=miss_val;
            indnull=logical(1-ind);
        else
            ind=1:m;
            indnull=[];
        end
        if Narg>3 %if using Binomial theory; rather, then empirical
            p=maf(i);
            mu=2*p;
            sig=sqrt(2*p*(1-p));
        else
            mu=mean(A(ind,i));
            sig=std(A(ind,i));
        end
        if sig~=0
            A(ind,i)=(A(ind,i)-mu)/sig;
        else
            A(ind,i)=0; %if variance is zero, then zero everything (due to order of limits)
        end
        A(indnull,i)=0; %zero missing
    end
    
    if Narg>1
        if dim==2 
            A=A';
        end
    end
end
function [B_final, F_final] = factorization_BF(X,A,k,W,options)

maxIter = 200;
if isfield(options, 'maxIter')
    maxIter = options.maxIter;
end
alpha = 1;
if isfield(options,'alpha')
    alpha = options.alpha;
end
%%%%%%===========Parameter settings===========
Norm = 2;
NormV = 1;


if min(min(A)) < 0
    error('Input should be nonnegative!');
end

[mFea,nSmp]=size(A);

%%%%%%=========== Weight matrix setting===========

if isfield(options,'weight') && strcmpi(options.weight,'NCW')
    feaSum = full(sum(X,2));
    D_half = (X'*feaSum).^.5;
    for i = 1:nSmp
        X(:,i) = X(:,i)/D_half(i);
    end
end

if isfield(options,'alpha_nSmp') && options.alpha_nSmp
    alpha = alpha*nSmp;    
end

W = alpha*W;% Weight matrix
DCol = full(sum(W,2));% Sum of row elements constitutes column vector DCol
D = spdiags(DCol,0,speye(size(W,1)));% Compose Diagonal D

%%%%%%%===========initialization================  
k1 = mFea;
if ~exist('U','var')
    B0 = abs(rand(k1,k));
    F0 = abs(rand(k,nSmp));
end

[B0,F0] = NormalizeUV(B0, F0', NormV, Norm);F0=F0';

Bk=B0;Fk=F0;
Akl=A;
Ak = A;
iter = 0; 
converged = 0;    
tol1=1e-5;tol2=1e-5;
%%%%%%%===========Update variables Z,A,B,F by iteration================

while ~converged  && iter < maxIter  
        iter = iter + 1;    
        alpha2=norm(Ak,'fro')/norm(Fk,'fro');
        %%%%%Update variables B and F, where Bk and Fk are the variables at the k-th iteration, 
        %%%%% and Bkl and Fkl are the variables at the k+1-th iteration.==================
        %%%========Update variables B==========
        Bkl=Bk.*(Akl*Fk')./((Bk*(Fk*Fk')));
        for i=1:size(Bkl,1)
            for j=1:size(Bkl,2)
                if Bkl(i,j)<0
                   Bkl(i,j)=0;
                else
                    Bkl(i,j)=Bkl(i,j);
                end
            end
        end
        %%%========Update variables F==========
        FF1=Bkl'*Akl+alpha2*Fk*W;
        FF2=(Bkl'*Bkl)*Fk+alpha2*Fk*D;
        Fkl=Fk.*(FF1)./(FF2);
        for i=1:size(Fkl,1)
            for j=1:size(Fkl,2)
                if Fkl(i,j)<0
                   Fkl(i,j)=0;
                else
                    Fkl(i,j)=Fkl(i,j);
                end
            end
        end
        
        [Bkl,Fkl] = NormalizeUV(Bkl, Fkl', NormV, Norm);
        Fkl=Fkl'; 
        Bwk = Bk;
        Fwk = Fk;
        Bk=Bkl;
        Fk=Fkl;
        %%%%%%%%%% Error
        Er1(iter,:)=abs(mean(A - Bk*Fk))./norm(Bk*Fk,'fro');
        er1(iter)=mean(Er1(iter,:)); 
        er=er1;
        temp = max ([norm(Bkl-Bwk,'fro'),norm(Fkl-Fwk,'fro')]);
        temp =temp/max(norm(X,'fro'));
        temp1 = norm( norm( (Ak - Bk*Fk),'fro'))/max(norm( Bk*Fk,'fro'));
        if temp1 < tol1 && temp < tol2
            converged = 1;
        end
end %end while
    
B_final = Bkl; %%% B_final  is finally B
F_final = Fkl; %%% F_final  is finally F

[B_final,F_final] = NormalizeUV(B_final, F_final', NormV, Norm);F_final=F_final';
t=1:iter;
figure
plot(t,er,'r-'),xlabel('Iteration times');ylabel('Error');
end



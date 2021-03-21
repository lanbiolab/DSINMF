function [A_final,B_final] = factorization_AB(X, k1, W, options)

maxIter = [];
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


if min(min(X)) < 0
    error('Input should be nonnegative!');
end

X=abs(X);
[n,m]=size(X);
A= abs(rand(n,k1));
B = abs(rand(k1,m));

[mFea,nSmp]=size(B);

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
L = D - W;
if isfield(options,'NormW') && options.NormW
    D_mhalf = DCol.^-.5;

    tmpD_mhalf = repmat(D_mhalf,1,nSmp);
    L = (tmpD_mhalf.*L).*tmpD_mhalf';
    clear D_mhalf tmpD_mhalf;

    L = max(L, L');
end

%%%%%%%===========initialization================


[A,B] = NormalizeUV(A, B', NormV, Norm);B=B';


Ak=A;Bk=B;
Ek=Bk;
Tk= zeros(k1,nSmp);
iter = 0; 
converged = 0;      
tol2=1e-5;
%%%%%%%===========add=================
O = zeros(n,m);%S矩阵的初始化
R = zeros(n,m); %用于非零位置的过滤矩阵
[x,y]=find(X == 0);
[k_,~]=size(x);
for i = 1:k_
   R(x(i),y(i)) = 1; 
end

% 计算S矩阵列向量的权重
col_sum = sum(X);
med = median(col_sum);
cof = col_sum/med;
P = diag(cof);

%%%%%%%===========Update variables A and B by iteration================

while ~converged  && iter < maxIter   

      iter = iter + 1;
      derta =5e+1;

      %%%%% Update variables A and B, where Ak and Bk are the variables at the k-th iteration, 
      %%%%% and Akl and Bkl are the variables at the k+1-th iteration.==================
      alpha1=norm(X+O,'fro')/norm(Ak,'fro');
      tau=norm(Ak*Bk-O,'fro')/norm(O,'fro');
      %%%========Update variables Z==========
      Ak1=Ak.*(((X+O)*Bk')./(Ak*(Bk*Bk')));
      %%%========Update variables A==========
      I=eye(k1);
      VV1=Ak1'*(X+O)+derta*Ek+Tk;
      VV2=((Ak1'*Ak1)+derta*I)*Bk;       
      Bkl=Bk.*(VV1./(VV2 ));
      % 更新S矩阵
      M=Ak1*Bkl-X;
      O=softx(M,2*tau);
      O=abs(O);
      O=O.*R;  
      % O=O*P;
      Ek=softx(Bkl-Tk/derta,2*alpha1/derta);
      Tk=Tk+1.618*derta*(Ek-Bkl);
    
      Awk = Ak;
      Bwk = Bk;
      Ak=Ak1;
      Bk=Bkl;
      %%%%%%%%%% Error
      Er1(iter,:)=abs(mean(X+O - Ak*Bk))./norm(Ak*Bk,'fro');
      er1(iter)=mean(Er1(iter,:));
      er=er1;
   
      temp =norm(X+O-Ak*Bk,'fro')/norm((X+O),'fro');
      if temp < tol2
        converged = 1;
      end

end %end while
  
A_final = Ak1; %%% Z_final  is finally Z
B_final = Bkl; %%% A_final  is finally A  
[A_final,B_final] = NormalizeUV(A_final, B_final', NormV, Norm);B_final=B_final';
t=1:iter;
figure(2);
plot(t,er,'r-'),xlabel('Iteration times');ylabel('Error');

end

function[y] = soft( x, T )
  if sum( abs(T(:)) )==0
       y = x;
  else
       y = max( abs(x) - T, 0);
       y = sign(x).*y;
   end
end    
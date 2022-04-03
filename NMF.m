function [W,H] = NMF(V,r)%设置分解矩阵的秩
[i,u]=size(V);  
W = rand(i,r);                            %初始化WH，为非负数
H = rand(r,u);
maviter=200;                                    %最大迭代次数
for iter=1:maviter
    W=W.*((V./(W*H))*H');           
    W=W./(ones(i,1)*sum(W));    
    H=H.*(W'*(V./(W*H)));
end

end
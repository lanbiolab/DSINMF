function [W,H] = NMF(V,r)%���÷ֽ�������
[i,u]=size(V);  
W = rand(i,r);                            %��ʼ��WH��Ϊ�Ǹ���
H = rand(r,u);
maviter=200;                                    %����������
for iter=1:maviter
    W=W.*((V./(W*H))*H');           
    W=W./(ones(i,1)*sum(W));    
    H=H.*(W'*(V./(W*H)));
end

end
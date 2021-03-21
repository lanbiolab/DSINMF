function [nmi,ami,ari,accuracy,glabel] = run_DSINMF(path_data,path_label,k_array_list)
    X=load(path_data);
    X = abs(X);
    X = X';%%%%%%%Rows are genes, columns are cell sample 
    X(all(X == 0,2),:)=[];
    real_label = load(path_label);

    %==============Constructing a weight matrix==============
    %Preset value before constructing weight matrix
    options = [];
    options.Metric = 'Cosine';
    options.NeighborMode = 'KNN';%KNN
    options.k =5;  %5 nearest neighbors
    options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 
    options.maxIter=400;
    
    options1 = [];
    options1.Metric = 'Cosine';
    options1.NeighborMode = 'KNN';%KNN
    options1.k =5;  %5 nearest neighbors
    options1.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 
    options1.maxIters=200;
    
    W = constructW(X',options);
    % W = SNN_Gauss(X', 10, 'cosine',1.0);
    X_in=X;
    if k_array_list < 2
        sprintf('k_value is not true!!!\n');
        exit(-1);
    end
    
    W = function_ds(W);
    [Z_final,A_final] = factorization_AB(X_in, k_array_list(1), W, options);
    A_in=A_final;
    for i=2:(length(k_array_list)+1)
        
        if i==(length(k_array_list)+1)
            ks = max(real_label);
        else
            ks = k_array_list(i);
        end
        [B,G_final] = factorization_BF(X_in,A_in,ks,W,options1);
        A_in=G_final;
    end
   

    l = zeros(1,size(G_final,2));
    for e=1:size(G_final,2)
        v=G_final(:,e);
        ma=max(v);
        [s,t]=find(v==ma);
        l(1,e)=s;
    end
    %%%%%%%%%%%%%%==================Performance evaluation===============================
    ll=real_label; %%%  the label originally identified by the authors
    l=l'; %%% Labels obtained by DRjCC
    [newl] = bestMap(ll,l); %% Permute label of l to match ll as good as possible
    glabel=newl;
    figure(1);
    mappedX = tsne(X');
    gscatter(mappedX(:,1), mappedX(:,2),newl);
    nmi=compute_NMI(ll,newl); %% Calculating the Normalized Mutual Information (NMI)
    ami=AMI(ll,newl); %% Calculating the Adjusted Mutual Information (AMI)
    ari = ARI(ll,max(ll),newl,max(newl)); %% Calculating the Adjusted Rand Index (ARI)
    pre_label =newl;
    if ~isempty(ll) 
    exact = find(pre_label == ll);
    accuracy = length(exact)/length(newl); %% Calculating the accuracy
    else
    accuracy = [];
    end
end

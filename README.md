
Detecting cell type from single cell RNA sequencing based on deep bi-stoachstic graph regularized matrix factorization
===========================================================================================================================
Overview:
----------------------------------------------------------------------------------------------------------------------------
###  This is code to do Detecting cell type from single cell RNA sequencing based on deep bi-stoachstic graph regularized matrix factorization given in the "experiment" section of the paper.The coding here is a generalization of the algorithm given in the paper. DSINMF is written in the MATLAB programming language. To use, please download the DSINMF folder and follow the instructions provided in the README.doc.
#### Files:
##### run_DSINMF.m - The main function.
##### factorization_AB.m - dimension reduction
##### factorization_BF.m - deep matrix factorization
##### constructW.m - Compute adjacent matrix W.
##### NormalizeUV.m - Normalize data.
##### bestMap.m - permute labels of L2 to match L1 as good as possible.
##### compute_NMI.m - Program for calculating the Normalized Mutual Information (NMI) between two clusterings.
##### AMI.m - Program for calculating the Adjusted Mutual Information (AMI) between two clusterings.
##### ARI.m - Program for calculating the Adjusted Rand Index ( Hubert & Arabie) between two clusterings.
##### hungarian.m - Solve the Assignment problem using the Hungarian method.
## Example:
### Follow the steps below to run DSINMF（also contained in the " excute_run_DSINMF.m" file）. Here use a real scRNA-seq data (Patel and Chu) set as an example.

```  
clc;
% dataSet={'Chu_ready','Patel_ready'};
dataSet={'Patel_ready'};  
k_array_list = [80,70,60];
disp(k_array_list);
c = 1; % 运行次数
cs=1; 
t1=clock;
for j = 1:c
    for i = 1:length(dataSet)           
            path_data = ['../DSINMF/',dataSet{i},'.txt'];
            path_label = ['../DSINMF/',dataSet{i},'_label.txt'];
            [nmi,ami,ari,accuracy] = run_DSINMF(path_data,path_label,k_array_list);
            %fprintf('data:%s\n',dataSet{i});
            fprintf('%f\t',nmi);
            fprintf('%f\t',ami);
            fprintf('%f\t',ari);
            fprintf('%f\t',accuracy);
    end
    fprintf('\n');  
 end 
t2=clock;
a=etime(t2,t1);
fprintf('%f\n',a); 
```

## Contact 
### Please send any questions or found bugs to lanwei@gxu.edu.cn 

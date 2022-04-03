clc;
% dataSet={'Chu_ready','Patel_ready'};
dataSet={'Darmanis'};
k_array_list = [80,70,60];
disp(k_array_list);
c = 10; % 运行次数
cs=1; 
t1=clock;

for j = 1:c
    for i = 1:length(dataSet)
            path_data = ['E:/single_data/featureSelection/',dataSet{i},'_data_p_value005.csv'];
            % path_data = ['E:/single_data/',dataSet{i},'_ready.txt'];
            path_label = ['E:/single_data/',dataSet{i},'_ready_label.txt'];
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





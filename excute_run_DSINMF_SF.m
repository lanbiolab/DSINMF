% clc;
% dataSet={'Chu_ready','Patel_ready'};
dataSet={'Haber'};
k_array_list = [80,70,60];
disp(k_array_list);
c = 1; % 运行次数
r = 1;
delta = 10;
for s = 1:r
    t1=clock;
    %selection_rate = (50-delta*s)/100;
    selection_rate = 0.7;
    % selection_rate = 1;
    for i = 1:length(dataSet)
        disp(dataSet{i});
        disp(selection_rate);
        for j = 1:c
            path_data = ['E:/single_data/',dataSet{i},'_ready.txt'];
            path_label = ['E:/single_data/',dataSet{i},'_ready_label.txt'];
            path_select = ['E:/single_data/featureSelection1/',dataSet{i},'_p_value_order.csv'];
            [nmi,ami,ari,accuracy] = run_DSINMF_SF(path_data,path_label,path_select,k_array_list,selection_rate);
            %fprintf('data:%s\n',dataSet{i});
            fprintf('%f\t',nmi);
            fprintf('%f\t',ami);
            fprintf('%f\t',ari);
            fprintf('%f\t',accuracy);
            fprintf('\n'); 
        end
    fprintf('\n');  
    end 
    t2=clock;
    a=etime(t2,t1);
    fprintf('%f\n',a);
end
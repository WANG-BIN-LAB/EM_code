clc,clear
filename='';
num = xlsread(filename,'Sheet1', 'A2:B361');
s=num(:,1);
systemByNode=num(:,2);
maindir=('E:\second_paper\result_all\high_low\P\RT_low');

for sub = 1:50 
for times = 1:50 
     FileName = ['modular_' num2str(times) 'subj_' num2str(sub) '.mat'];
     fulldir= fullfile(maindir,'\',FileName);
    a = load(fulldir);
    MA1 = a.mat;
    for i = 1:length(s)
        for j = 1:length(s)
            MA(i,j) = MA1(s(i),s(j));
        end
    end
    I = integration(MA,systemByNode);

    for_I_mean(times,:) = I';
    R = recruitment(MA,systemByNode);
    for_R_mean(times,:) = R';
end
    mean_times_network_integration(sub,:) = mean(for_I_mean);
    mean_times_network_recruitment(sub,:) = mean(for_R_mean);
end
xlswrite('',mean_times_network_integration);
xlswrite('',mean_times_network_recruitment);








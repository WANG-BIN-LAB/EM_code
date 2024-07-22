clc;clear;close all
N=360;%N_nodes
L=228;%time_points
maindir='';  
for sub = 1:50
    for times = 1:50 
        FileName = ['S_' num2str(times) 'th_times_' num2str(sub) 'subject.mat'];
        fulldir= fullfile(maindir,'\',FileName);
        S=load(fulldir);
        A=S.S;
        [mat] = community(A,N,L);
        savefile='';    
        save(strcat(savefile,'modular_',num2str(times),'subj_',num2str(sub),'.mat'),'mat'); 
    end
end

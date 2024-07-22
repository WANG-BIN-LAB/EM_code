
clc,clear
% The folder path where the data is located
folder=''; 
%Find all mat files
files=dir([folder,'\*.mat']); 
%Initialize the cellular array and store the original matrix sequence
ROISignals_1=cell(1,1); 
for i=1:81 
ROISignals_1{i}=cell2mat(struct2cell(load(strcat('',files(i).name))));
end
for m=1:81
fMRI_data_all_subjects(m,:,:)=ROISignals_1{1,m};
end
%save
save('\fMRI_data_all_subjects.mat','fMRI_data_all_subjects')
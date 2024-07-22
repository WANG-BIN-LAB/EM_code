%%
clc,clear
%%
F=xlsread('pamret_scores_of_mean_integration.xlsx','sheet1','B1:B81');
M=xlsread('pamret_scores_of_mean_integration.xlsx','sheet1','C1:MX81');
hcp_360=xlsread('E:\second_paper\result_all\memory_scores\360_7.xlsx','sheet1','A1:B360');
%%
beishi=81;

divide_result_M=zeros(7,beishi,360);
P=zeros(7,1);
%mu=[6,54,39,56,23,23,50,15,77,7,4,6];
mu=[59;53;44;48;29;45;82];
for i=1:7
    for j=1:360
        if hcp_360(j,2)==i
            divide_result_M(i,:,j)=M(:,hcp_360(j,1));
        end
    end
    
    %pamenc_integration_F=zeros(1,beishi);
    pamenc_integration_M=zeros(1,beishi);
    
%     for Z=1:beishi
%         pamenc_integration_F(1,Z)=sum(divide_result_F(i,Z,:))/mu(i);
%     end
    for Z=1:beishi
        pamenc_integration_M(1,Z)=sum(divide_result_M(i,Z,:))/mu(i);
    end
    %
    p=corr(F,pamenc_integration_M');
    P(i,1)=p;
end
%%
%BRAIN
pamenc_integration_zong=zeros(beishi,1);
for i=1:beishi
    pamenc_integration_zong(i,1)=mean(M(i,:));
end
p_quan=corr(F,pamenc_integration_zong);



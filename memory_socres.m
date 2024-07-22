clc,clear
file_path = '';
files = dir(fullfile(file_path, '*.tsv')); %��������tsv�ļ�

%trial_type CONTROL��Ϊ����
type_1='CORRECTLY';
type_2='INCORRECTLY';
%response_str
response_str_1='SURE_CORRECT';%4 �� %1
response_str_2='MAYBE_CORRECT';%3 �� %2
response_str_3='MAYBE_INCORRECT';%2 �� %3
response_str_4='SURE_INCORRECT';%1 �� %4 
result=zeros(81,1);
for i=1:81
    all=zeros(104,1);
    data = readtable(strcat('E:\second_paper\result_all\memory_scores\task_pamret_events\',files(i).name), 'FileType', 'text', 'Delimiter', '\t');
    % 'CORRECTLY'�ж�
    for j=1:104
        if isequal(cell2mat(data(j,4).trial_type(1,:)),type_1)
            if isequal(cell2mat(data(j,7).response_str(1,:)),response_str_1)
                all(j,1)=4;
            elseif isequal(cell2mat(data(j,7).response_str(1,:)),response_str_2)    
                all(j,1)=3;
            elseif isequal(cell2mat(data(j,7).response_str(1,:)),response_str_3)   
                all(j,1)=2;
            elseif isequal(cell2mat(data(j,7).response_str(1,:)),response_str_4)   
                all(j,1)=1;
            end  
         % 'INCORRECTLY'�ж�
        elseif isequal(cell2mat(data(j,4).trial_type(1,:)),type_2)
             if isequal(cell2mat(data(j,7).response_str(1,:)),response_str_1)
                all(j,1)=1;
            elseif isequal(cell2mat(data(j,7).response_str(1,:)),response_str_2)    
                all(j,1)=2;
            elseif isequal(cell2mat(data(j,7).response_str(1,:)),response_str_3)   
                all(j,1)=3;
            elseif isequal(cell2mat(data(j,7).response_str(1,:)),response_str_4)   
                all(j,1)=4;
            end             
        end
    end
    result(i,1)=sum(all(:,1))/80;
    
end



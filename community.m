function [ mat ] = community( A,N,L ) 

mat=zeros(N,N);
%for l=1:L
for k=1:L   
    for i=1:N
        for j=1:N
            if A(i,k)==A(j,k)
                mat(i,j)=mat(i,j)+1;
            else
                mat(i,j)=mat(i,j)+0;
            end
        end
    end  
end
 mat=mat/L;

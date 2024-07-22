function S=postprocess_ordinal_multilayer(S,T,max_coms,verbose)
% POSTPROCESS_ORDINAL_MULTILAYER post-process an ordered multilayer partition 
% to maximise persistence when using uniform ordinal coupling with the 
% multilayer quality function in Mucha et al. 2010. 
%在Mucha等人2010年使用与多层质量函数一致的有序耦合时，对有序的多层分区进行后处理以最大化持久性。
% Version: 2.1.1
% Date: Mon 27 Feb 2017 19:15:16 EST
% 
% (see Bazzi et al. 2016 for more detail on persistence and
% post-processing when using the multilayer quality function in Mucha et
% al. 2010 for ordered layers).
%
% Call as:
%
%     S=postprocess_ordinal_multilayer(S)
%
%     S=postprocess_ordinal_multilayer(S,T)
%
%     S=postprocess_ordinal_multilayer(S,T,max_coms)
%
% Input: 
% 
%     S: multilayer partition
%
%     T: number of layers of S (defaults to 'size(S,2)')
%
%     max_coms: only run function when input partition has less than 
%         'max_coms' communities, otherwise return input partition.
%         (defauts to 'inf')
%
% Output:
%
%     S: post-processed multilayer partition 
%
% Post-processes a multilayer partition to maximise persistence without 
% changing the community structure within layers, using the Hungarian
% algorithm to solve the optimal assignment problem for each consecutive
% pair of layers. 
%对一个多层分区进行后处理，在不改变层内社区结构的情况下最大化持久性，使用Hungarian算法解决每对连续层的最优分配问题。
%Note that this procedure always increases
% multilayer modularity with ordinal uniform interlayer coupling for any 
% non-zereo value of coupling strength omega. This function can be 
% particularly useful when using the multilayer quality function in 
% Mucha et al. 2010 with low values of ordinal uniform interlayer coupling. 
%注意，对于耦合强度omega的任何非零值，这个过程总是通过有序均匀层间耦合来增加多层模块化。这个函数在Mucha等人2010年使用低阶均匀层间耦合值的多层质量函数时特别有用。
%   References:
%
%     Mucha, Peter J., Thomas Richardson, Kevin Macon, Mason A. Porter, and
%     Jukka-Pekka Onnela. "Community Structure in Time-Dependent,
%     Multiscale, and Multiplex Networks," Science 328, 876-878 (2010).
%
%     Bazzi, Marya, Mason A. Porter, Stacy Williams, Mark McDonald, Daniel
%     J. Fenn, and Sam D. Howison. "Community Detection in Temporal 
%     Multilayer Networks, with an Application to Correlation Networks", 
%     MMS: A SIAM Interdisciplinary Journal 14, 1-41 (2016).
%
%   Citation: If you use this code, please cite as
%       Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla and Peter J. Mucha,
%       "A generalized Louvain method for community detection implemented in
%       MATLAB," http://netwiki.amath.unc.edu/GenLouvain (2016).

if nargin<2||isempty(T)
    T=size(S,2);
end



N=numel(S)/T;


if nargin<3||isempty(max_coms)
    max_coms=inf;
end

if nargin<4||isempty(verbose)
    verbose=false;
end
if verbose
    mydisp=@(s) disp(s);
else
    mydisp=@(s) [];
end

[N0,T0]=size(S);


if max(S(:))<max_coms % don't do anything if too many communities for performance
    S=reshape(S,N,T); 
    if verbose
        p0=ordinal_persistence(S);
    end
    max_com=max(S(:,1));
    for i=2:T
        [u1,~,e1]=unique(S(:,i-1)); % unique communities in previous layer
        [u2,~,e2]=unique(S(:,i)); % unique communities in this layer
        G1=sparse(e1,1:length(e1),1); % community assignment matrix for previous layer 前一层的社区分配矩阵
        G2=sparse(1:length(e2),e2,true); % community assignment matrix for this layer 这一层的社区分配矩阵
        overlap=G1*G2; % node overlap matrix between communities in the two layers 两层中社区之间的节点重叠矩阵
        dist=sum(overlap(:))-overlap;
        S2=assignmentoptimal(dist'); % find best assignment for communities in current layer 为当前层的社区寻找最佳分配

        for j=1:length(u2)
            if S2(j)~=0&&overlap(S2(j),j)
                S(G2(:,j),i)=u1(S2(j)); % update assignment
            else
                S(G2(:,j),i)=max_com+1; % get new label if not assigned
                max_com=max_com+1;
            end
        end
    end
    if verbose
        p1=ordinal_persistence(S);
        fprintf('Improvement in persistence: %g\n',p1-p0);
    end
else
    mydisp('number of communities exceeds ''max_coms'', skipping postprocessing')
end

% return in original format
S=reshape(S,N0,T0);
end


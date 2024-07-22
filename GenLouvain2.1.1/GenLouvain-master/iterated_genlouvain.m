function [S,Q,n_it]=iterated_genlouvain(B,limit,verbose,randord,randmove,S0,postprocessor)
% Optimise modularity-like quality function by iterating GenLouvain until convergence.
% (i.e., until output partition does not change between two successive iterations)
% 通过迭代GenLouvain直到收敛来优化模块化的质量函数。
%(即。，直到两个连续迭代之间的输出分区不变)
% Version: 2.1.1
% Date: Mon 27 Feb 2017 19:15:16 EST
%
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B) with matrix B repeatedly implements 
%   GenLouvain until convergence to an output partition.矩阵B反复实现GenLouvain，直到收敛到一个输出分区。
%   After the first 
%   iteration, the starting partition of each subsequent iteration is the 
%   output partition of the previous iteration. The output vector S encodes 
%   the community assignments of the output partition after the final 
%   iteration, with S(i) identifying the community to which node i has been 
%   assigned.
%在第一次迭代之后，每个后续迭代的起始分区是前一次迭代的输出分区。
%输出向量S在最后一次迭代后编码输出分区的社区分配，S(i)标识节点i被分配给哪个社区。
%The output Q gives the quality of the partition S of the 
%   network. The output n_it gives the observed number of iterations until 
%   convergence.
%输出Q给出了网络分区的质量。输出n_it给出了在收敛之前观察到的迭代次数。
%   NOTE: The matrix represented by B must be both symmetric and square. 
%   This condition is not checked thoroughly if B is a function handle, but
%   is essential to the proper use of this routine. When B is a matrix, 
%   non-symmetric input is symmetrised (B=(B+B')/2), which preserves the
%   quality function.
%由B表示的矩阵必须是对称的方阵。如果B是一个函数句柄，则不会彻底检查这个条件，但是对于正确使用这个例程是必要的。当B为矩阵时，将非对称输入将对称处理为(B=(B+B')/2)，保留质量函数。
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B) with function handle B such that 
%   B(i) returns the ith column of the modularity/quality matrix uses this
%   function handle (to reduce the memory footprint for large networks) 
%   until the number of groups is less than 10000 and then builds the B 
%   matrix corresponding to the new aggregated network in subsequent passes.  
%   Use [S,Q,n_it] = ITERATED_GENLOUVAIN(B,limit) to change this 
%   default=10000 limit.
%
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B,limit,0) suppresses displayed text 
%   output.
%
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B,limit,verbose,0) forces index-ordered 
%   (cf. randperm-ordered) consideration of nodes, for deterministic results 
%   with randord = 'move'.
%
%   [S,Q,n_it]=ITERATED_GENLOUVAIN(B,limit,verbose,randord,randmove) controls 
%   additional randomization to obtain a broader sample of the quality 
%   function landscape. The possible values for 'randmove' are
%       'move': always move node under consideration to the community that 
%           results in maximal improvement in modularity (default)
%总是考虑将节点移动到社区，以最大限度地提高模块化(默认)
%       'moverand': move the node under consideration to a community chosen
%           uniformly at random from all moves that increase the qualilty 
%           function
%将考虑中的节点从所有增加质量函数的移动中随机均匀地移动到一个社区
%       'moverandw': move the node under consideration to a community chosen
%           at random from all moves that increase the qualilty where the
%           probability of choosing a particular move is proportional to
%           its increase in the quality function
%将考虑中的节点移动到随机选择的社区中，从所有增加质量的移动中选择一个特定移动的概率与质量函数的增加成正比
%       0: equivalent to 'move' (provided for backwards compatibility)
%       1: equivalent to 'moverand' (provided for backwards compatibility)
%
%   'moverand', and 'moverandw' mitigate some undesirable behavior for 
%   "multilayer" modularity with ordinal coupling ('moverandw' tends to be 
%   better behaved for large values of the interlayer coupling). 
%'moverand'和'moverandw'减轻了使用有序耦合的“多层”模块化的一些不良行为('moverandw'对于较大的层间耦合值表现得更好）
%   With 'move', the algorithm exhibits an abrupt change in behavior when the 
%   strength of the interlayer coupling approaches the maximum value of the
%   intralayer modularity matrices (see Bazzi et al. 2016 for more detail).
%对于“move”，当层间耦合的强度接近层内模块化矩阵的最大值时，该算法的行为会发生突变(详见Bazzi等人2016年的文章)。
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B,limit,verbose,randord,randmove,S0) 
%   uses S0 as an inital partition for the first iteration (the starting 
%   partition of subsequent iterations is the output partition of the previous
%   iteration). 
%使用S0作为第一次迭代的初始分区(后续迭代的起始分区是前一次迭代的输出分区)。
%The default choice for S0 is all singletons (and given by a 
%   length(B) by 1 vector). If a user specifies S0, it needs to satisfy 
%   numel(S0) = length(B). Note that the size of S will match the size of S0. 
%   In a multilayer setting, a user may have to reshape S appropriately 
%   (e.g., reshape(S,N,T), where N is the number of nodes in each layer and 
%   T is the number of layers).
%S0的默认选择是所有的单个节点(由length(B)*1向量给出)。如果用户指定了S0，它需要满足numel(S0) = length(B)。（S0的总元素个数=B的列数）
%注意，S的大小将与S0的大小匹配。在多层设置中，用户可能需要对S进行适当的reshape(例如，reshape(S,N,T)，其中N是每层的节点数，T是层数)。
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B,limit,verbose,randord,randmove,S0,
%   postprocessor) applies the function handle postprocessor to the
%   output partition of GENLOUVAIN at each iteration. The function handle
%   postprocessor takes a partition as input and returns a partition as 
%   output. 在每次迭代时将函数句柄postprocessor应用于GENLOUVAIN的输出分区。函数句柄 postprocessor用一个分区作为输入，并返回一个分区作为输出。
%To ensure convergence, the function postprocessor should
%   increase the value of the (multilayer) quality function for a partition 
%   (e.g, kernigan-lin type algorithms in a monolayer or multilayer setting, 
%   and "postprocess-ordinal-multilayer.m",
%   "postprocess-categorical-multilayer.m" in HelperFunctions for a multilayer 
%   setting) 
%为了保证收敛，函数postprocessor应该为一个分区增加(多层)质量函数的值(例如在单层或多层设置以及 在HelperFunctions的多层设置中"postprocess-ordinal-multilayer.m",
%   "postprocess-categorical-multilayer.m"的kernigan-lin类型算法)
%   Example on multilayer network quality function of Mucha et al. 2010
%   (using multilayer cell A with A{s} the adjacency matrix of layer s)
%
%   gamma = 1; omega = 0.1; 
%   N=length(A{1});
%   T=length(A);
%
%   B = multiord(A,gamma,omega); % multiord.m computes multilayer
%   modularity matrix B with homogeneous ordinal interlayer coupling w and 
%   a Newman-Girvan null model on each layer (see more detail in
%   documentation of multiord.m)
%
%   PP = @(S) postprocess_ordinal_multilayer(S,T); % define postprocessing
%   function handle that increases multilayer modularity without changing 
%   intralayer partitions in an ordered multilayer networks (see 
%   postprocess_ordinal_multilayer.m for more detail)
%定义后处理函数句柄，在不改变有序多层网络的层内分区的情况下增加多层模块化(参见postprocess_ordinal_multilayer)。m表示更多细节)
%   [S,Q,n_it] = ITERATED_GENLOUVAIN(B,10000,0,1,'moverandw',[], PP);
%   S = reshape(S, N, T);
%
%   Notes:
%
%     Under default options, this routine can return different results from
%     run to run because it considers nodes in pseudorandom (randperm)
%     order.  Because of the potentially large number of nearly-optimal
%     partitions (Good et al. 2010), one is encouraged to investigate
%     results of repeated applications of this code (and, if possible, of
%     other computational heuristics).  To force deterministic behavior with 
%     randord = 'move', ordering nodes by their index, pass zero as the 
%     fourth input: ITERATED_GENLOUVAIN(B,limit,verbose,0).
%在默认选项下，这个例程可以在每次运行时返回不同的结果，因为它以伪随机(randperm)顺序考虑节点。由于可能存在大量接近最优的分区(Good et al. 2010)，因此鼓励研究这段代码重复应用的结果(如果可能，还可以研究其他计算启发法)。
%要使用randord = 'move'强制执行确定性行为，按节点的索引对节点排序，传递0作为第四个输入:ITERATED_GENLOUVAIN(B,limit,verbose,0)。
%     This algorithm is only "Louvain-like" in the sense that the two
%     phases are used iteratively in the same manner as in the Louvain
%     algorithm (Blondel et al. 2008).  Because it operates on a general
%     quality/modularity matrix B, it does not include any analytical
%     formulas for quickly identifying the change in modularity from a
%     proposed move nor any improved efficiency obtained by their use.  If
%     your problem uses one of the well-used null models included in other
%     codes, those codes should be much faster for your task.
%该算法只是在两个阶段的迭代使用方式与Louvain算法相同的情况下才具有“Louvain样”(Blondel et al. 2008)。
%因为它是在一个通用的质量/模块化矩阵B上运行的，所以它不包括任何快速识别模块化变化的分析公式，也不包括通过使用它们而获得的任何效率的提高。
%如果您的问题使用了其他代码中包含的常用空模型，那么这些代码对于您的任务来说应该会快得多。
%     Past versions had a problem where accumulated subtraction error might
%     lead to an infinite loop with each pass oscillating between two or
%     more partitions yet incorrectly identifying increases in quality.  We
%     believe this problem has been corrected by the relative change checks
%     in lines 178 and 269.  If you encounter a similar problem, notify
%     Peter Mucha (<a href="mailto:mucha@unc.edu">mucha@unc.edu</a>).
%
%     The output Q provides the sum over the appropriate elements of B
%     without any rescaling.  输出Q提供了B的适当元素的和，而不需要任何重新缩放。
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
%
%   References:
%     Blondel, Vincent D., Jean-Loup Guillaume, Renaud Lambiotte, and
%     Etienne Lefebvre, "Fast unfolding of communities in large networks,"
%     Journal of Statistical Mechanics: Theory and Experiment, P10008
%     (2008).
%
%     Fortunato, Santo, "Community detection in graphs," Physics Reports
%     486, 75-174 (2010).
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
%     Porter, M. A., J. P. Onnela, and P. J. Mucha, "Communities in
%     networks," Notices of the American Mathematical Society 56, 1082-1097
%     & 1164-1166 (2009).
%
%   Acknowledgments:
%     A special thank you to Stephen Reid, whose greedy.m code was the
%     original version that has over time developed into the present code. 
%
%     Thank you also to Dani Bassett, Jesse Blocher, Mason Porter and Simi
%     Wang for inspiring improvements to the code.
%
%   Citation: If you use this code, please cite as
%       Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla and Peter J. Mucha,
%       "A generalized Louvain method for community detection implemented in
%       MATLAB," http://netwiki.amath.unc.edu/GenLouvain (2016).
%
%   See also genlouvain HelperFunctions


%set default for maximum size of modularity matrix
if nargin<2
    limit = [];
end

%set level of reported/displayed text output
if nargin<3||isempty(verbose)%判断verbose是否为空，如果为空，结果为1，否则为0.
    verbose = false;
end


%set randperm- v. index-ordered
if nargin<4
    randord = [];
end

%set move function (maximal (original Louvain) or random improvement)
if nargin<5
    randmove=[];
end

% set initial partition
if nargin<6
    S0=[];
end

% set postprocessing function
if nargin<7
    postprocessor=[];
end

% verbose output switch
if verbose
    mydisp = @(s) disp(s);
else
    mydisp = @(s) [];
end

S_old=[];
n_it=1;
mydisp('Iteration 1');
[S,Q]=genlouvain(B,limit,verbose,randord,randmove,S0);

mydisp('');

Q_old=-inf;
while ~isequal(S,S_old)&&(Q-Q_old)>10*eps
    n_it=n_it+1;
    S_old=S;
    Q_old=Q;
    
    mydisp(sprintf('Iteration %u',n_it));
    if ~isempty(postprocessor)
        S=postprocessor(S);
    end
    [S,Q]=genlouvain(B,limit,verbose,randord,randmove,S);
    mydisp(sprintf('Improvement in modularity: %f\n',Q-Q_old));
end

end

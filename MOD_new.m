%% preliminaries
clc,clear
%load 3D fMRI data; in our analysis the size of the 3D array was [1003,4800,25] (subjects,time-points,nodes)
load(''); 
n_subjects = size(fMRI_data_all_subjects,1); % total number of subjects
time_points = size(fMRI_data_all_subjects,2); % total number of time-points in fMRI data
N_nodes = size(fMRI_data_all_subjects,3); % number of brain nodes in fMRI data
window_overlap = 1; % overlap between sliding-windows in number number of TRs
window_length =15; % total window length (here, 139 time-points correspond to 100s with a TR of 0.72s)
intralayer_resolution = 1; % intralayer resolution for multilayer modularity
interlayer_resolution = 1; % interlayer resolution for multilayer modularity

%% pre-allocate arrays - this is the output data (also saved at the end of script)
network_switching = zeros(n_subjects,N_nodes); % network switching data
Q_value = zeros(1,n_subjects); % Q-value of networks
number_of_communities = zeros(n_subjects,floor((time_points-window_length)/window_overlap)+1); % number of communities
modularity_iterations = zeros(1,n_subjects); % number of iterations needed before the multilayer modularity algorithm convergestest1

%  parallel_loop_info = parpool(); % this will intitialize the parallel pool used next
%  parallel_loop_NumWorkers = parallel_loop_info.NumWorkers; % number of workers are needed to estimate total time-needed to run this script (see line #202)

for times=1:50
for n = 1:n_subjects
    %% Time-varying connectivity (A is for dynamic connectivity; AA is for multilayer modularity)
    A = zeros(floor((time_points-window_length)/window_overlap)+1,N_nodes,N_nodes); % pre-allocate connectivity matrices - multilayer modularity
    AA = cell(1,floor((time_points-window_length)/window_overlap)+1); % pre-allocate connectivity matrices - dynamic connectivity
    fMRI_data_current_subject = squeeze(fMRI_data_all_subjects(n,:,:)); % extract fMRI data from subject n
    
    for i = 1:floor((time_points-window_length)/window_overlap)+1
        current_window = (squeeze(fMRI_data_current_subject((i:i+window_length-1)*window_overlap,:))).*repmat(hamming(window_length),1,N_nodes); % current fMRI window - use hamming window
        %current_window = squeeze(fMRI_data_current_subject((i:i+window_length-1)*window_overlap,:)); % current fMRI window - use hamming window

        sliding_window = corrcoef(current_window); % correlation-based sliding window analysisw on fMRI data
     
        A(i,:,:) = sliding_window; % matrices for dynamic connectivity analysis
        sliding_window(sliding_window<0) = 0; % remove negative correlation values for multilayer modularity analysis

        AA{i} = sliding_window; % matrices for multilayer modularity analysis
      
    end
    
    %% Ordinal iterative multilayer modularity calculation
    [B,mm] = multiord(AA,intralayer_resolution,interlayer_resolution); % B=supra-adjacency matrix
    PP = @(S) postprocess_ordinal_multilayer(S,floor((time_points-window_length)/window_overlap)+1); % multilayer modularity set-up 
    [S,Q1,mod_iter_tmp] = iterated_genlouvain(B,10000,0,1,'moverandw',[],PP); % S = modularity assignments; Q1 = modularity value; mod_iter_tmp = iterations
    S = reshape(S,N_nodes,floor((time_points-window_length)/window_overlap)+1); % 2D representation of modularity assignments
    
    savefile='';
    save(strcat(savefile,'S_',num2str(times), 'th_times_',num2str(n),'subject.mat'),'S'); 
    
%   
 end
% %% Save output data
% save(['switchings_analysis with_' num2str(times) '_th_' num2str(N_nodes) '_nodes_' num2str(time_points) '_timepoints_and_' num2str(n_subjects) '_subjects.mat']);
 end
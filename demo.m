clear all
addpath(genpath('libs'))

%% Load shape
load('data\cat0_rem5000.mat')
M1 = M;
load('data\cat5_rem5000.mat')
M2 = M;

%% compute spectral decomposition of the shape laplacian
n_eigenvalues = 100;
[M1.evecs, M1.evals]=calc_eig_FEM(M1,n_eigenvalues);
[M2.evecs, M2.evals]=calc_eig_FEM(M2,n_eigenvalues);


%% comput MMS point signature
tic
M1.MMS = calc_MMS(M1.evecs, M1.evals);
M2.MMS = calc_MMS(M2.evecs, M2.evals);
toc

%% nearest neighbor match and error
nm = 500;
match = knnsearch(M1.MMS,M2.MMS,'k',nm);
gt_match = M1.remesh_idx(M2.orig_idx)';

hit = match == gt_match;
perc_match = 100*(1:nm)/M1.n;
hit_rate = [];
for i=1:nm
    hit_rate(i) = mean(any(hit(:,1:i),2));
end
plot(perc_match, hit_rate);
xlabel '% matches'
ylabel 'hit ratio'
grid on






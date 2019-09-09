function [Tevecs,Tevals,Sl] = calc_eig_FEM(T,k)
%%
[W,Sc,Sl] = calc_LB_FEM(T);
Si = sparse(diag(sqrt(1./diag(Sl))));
LAP = (Si*W*Si);
LAP = (LAP+LAP')./2;

[Tevecs_,Tevals_] = eigs(LAP,k,-1e-5);
[Tevecs,Tevals] = eigs(W,Sl,k,-1e-5);
Tevals = diag(Tevals);

[Tevals, idx] = sort(Tevals);
Tevecs = Tevecs(:,idx);


% Tevals2 = eig(LAP); 
% Tevals2 = Tevals2(1:k);
% close all
% plot(Tevals)
% hold on 
% plot(Tevals2)

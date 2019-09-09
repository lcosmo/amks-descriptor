function [MMSdiag] = calc_MMS(evec,eval,time, desc_size, variance)



%% set default parameters
if nargin<3; time = 0.5; end
if nargin<3; desc_size = 100; end
if nargin<3; variance = 6; end

n = size(evec,1);

%% check if GPU is available
try
   gpuArray(1);
   canUseGPU=true;
catch
 canUseGPU=false;
end

%% Compute bands
log_E = log(max(eval, 1e-6))';
e = linspace(log_E(2), (max(log_E))/1.02, desc_size);
sigma = (e(2)-e(1))*variance;

%% Compute descriptor
P=evec.^2;
MMSdiag = zeros(size(evec,1),desc_size);

if canUseGPU 
    P = gpuArray(P); 
    MMSdiag = gpuArray(MMSdiag);
end

F1 = sinc(time*(eval-eval')./pi);
for j=1:desc_size
    %%
    F2 = exp((-(e(j) - log_E).^2) ./ (2*sigma.^2));
    F2 = F2.*F2';
    F = F1.*F2./sum(triu(F2),1:2);
    MMSdiag(:,j) = dot(P',(P*F)');
end

if canUseGPU 
    MMSdiag = gather(MMSdiag);
end

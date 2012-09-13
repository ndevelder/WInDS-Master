clear all

addpath(genpath(fullfile(cd))); %Add directories to search path

% load('04-06-2012_1127_NRELrotor_static_rated.mat', '-mat');

load('05-03-2012_1702_NRELrotor_static_rated.mat', '-mat');


p = 100;
time = 5;

P = wake.domain{p};
F1 = wake.domain{p}(:,:,2:end,:);
F2 = wake.domain{p}(:,:,1:end-1,:);
gamma = wake.gamma.trail{p};
rc = wake.rc_eff.trail{p};


% P = pos.bound(:,:,time,:);
% F1 = wake.domain{time}(1:end-1,:,:,:);
% F2 = wake.domain{time}(2:end,:,:,:);
% gamma = wake.gamma.shed{time};
% rc = wake.rc_eff.shed{time};

gpu = [128 1 0];
gpuuse = 'true';
co = 1000.0;
casetype = 'fuln';
cmod = 'none';
d = 0.01;
uind = zeros(size(P));

tic;
[uind L] = BiotSavartMex(F1,F2,P,gamma,rc,d,cmod,co,casetype,gpuuse, gpu);
mex = toc;

% tic;
% [trueuind trueL] = BiotSavart(F1,F2,P,gamma,rc,d,co,'full');
% mat = toc;
% 
% differ = sum(sum(sum(sum(trueuind-uind))));
% 
% a = sprintf('Mex time: %3.9f s Matlab time: %3.9f s Summed Difference: %3.9f\n',mex,mat,differ);
% 
% disp(a)

%Uncomment when using NVPP Profiler
exit
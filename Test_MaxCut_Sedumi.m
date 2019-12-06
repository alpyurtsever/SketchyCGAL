%%  Test setup for MaxCut SDP - Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)
% Solved with Sedumi
% To replicate the results in Figure 7.1 in the paper

%% Choose and load data
% NOTE: You need to download data from GSET and locate them to under the
% "FilesMaxCut/data/G/" folder (resp. DIMACS10, "FilesMaxCut/data/DIMACS10/").

maxcut_data = 'G/G1'; % you can choose other data files as well
% maxcut_data = 'DIMACS10/belgium_osm';

%% Preamble
rng(0,'twister');

% Add SEDUMI to path!
% NOTE: You need to download SEDUMI if you don't already have it.

% addpath(genpath('DESTINATION OF SEDUMI'));

%% load data

% load(['./FilesMaxCut/data/',maxcut_data]);
load(['/home/yurtseve/MATLAB_COMMON/DATASETS/MaxCut/',maxcut_data]);

n = size(Problem.A,1);
L = spdiags(Problem.A*ones(n,1),0,n,n) - Problem.A;
L = 0.5*(L+L');
clearvars Problem;

%% Construct Sedumi Problem

K.f = 0;
K.l = 0;
K.q = [];
K.r = [];
K.s = n;
K.scomplex = [];
K.ycomplex = [];
c = (-0.25).*L(:);
clearvars L;
indI = (0:n-1)'*(n+1) + 1;
indJ = (1:n)';
valA = ones(n,1);
b = ones(n,1);
At = sparse(indI,indJ,valA,n^2,n);
clearvars indI indJ valA;

% pars.eps = 1e-4;

%% Start memory logging
% NOTE: This works only on Unix systems. 

hmL = memLog([maxcut_data,'sedumi']);
hmL.start;
MEMBEGIN = hmL.prompt;

%% Solve with SEDUMI

timer = tic;
cputimeBegin = cputime;

[x,y,info] = sedumi(At,b,c,K);
% [x,y,info] = sedumi(At,b,c,K,pars);

cputimeEnd = cputime;
totalTime = toc(timer);

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

%% Stop memory logging

MEMEND = hmL.prompt;
hmL.stop;
out.memory = (MEMEND - MEMBEGIN)/1000; %in MB

%% Compute cutvalue

R = 10; % rank/sketch size parameter

cutvalue = 0;
c = reshape(c,[n,n]);
x = reshape(x,[n,n]);
x = 0.5*(x+x'); % sedumi does not return symmetric outputs
[u,~] = eigs(x, R, 'LM');
for t = 1:R
    sign_evec = sign(u(:,t));
    rankvalue = -(sign_evec'*(c*sign_evec));
    cutvalue = max(cutvalue, rankvalue);
end

out.cutvalue = cutvalue;
out.primalObj = c(:)'*x(:);
out.primalFeas = norm(diag(x)-b)/max(norm(b),1);
out.sdpFeas = abs(min(min(eig(x)),0));

%% Save Results

if ~exist(['results/MaxCut/',maxcut_data],'dir'), mkdir(['results/MaxCut/',maxcut_data]); end
save(['results/MaxCut/',maxcut_data,'/SeDuMi.mat'],'info','out','-v7.3');

%% Last edit: Alp Yurtsever - November 18, 2019
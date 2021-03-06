%%  Test setup for MaxCut SDP - Solved with MoSeK
%% Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)

%% Choose data
% NOTE: You need to download data from GSET and locate them to under the
% "FilesMaxCut/data/G/" folder (resp. DIMACS10, "FilesMaxCut/data/DIMACS10/").

% maxcut_data = 'G/G1'; % you can choose other data files as well
% maxcut_data = 'DIMACS10/delaunay_n10';

%% Preamble
rng(0,'twister');

% Add Mosek to path!
% NOTE: You need to download MoSeK if you don't already have it.
% MoSeK requires a license, personal academic license for MoSeK is free.

% addpath(genpath('DESTINATION OF MOSEK'));

%% Load data

load(['./FilesMaxCut/data/',maxcut_data]);

n = size(Problem.A,1);
C = spdiags(Problem.A*ones(n,1),0,n,n) - Problem.A;
clearvars Problem;
C = 0.5*(C+C'); % symmetrize if not symmetric
C = (-0.25).*C;
b = ones(n,1);

%% Construct MoSeK problem

% Loss
prob.bardim = size(C,1);
[prob.barc.subk, prob.barc.subl, prob.barc.val] = find(tril(C));
prob.barc.subj = ones(size(prob.barc.subk));
 
% diag(X) = 1 constraint
prob.a 			= sparse(n,n);
prob.bara.subi 	= (1:n)'; 
prob.bara.subj 	= ones(n,1);
prob.bara.subk 	= (1:n)'; 
prob.bara.subl 	= (1:n)'; 
prob.bara.val  	= ones(size(prob.bara.subi)); 
prob.blc 		= ones(size(prob.bara.subi)); 
prob.buc 		= ones(size(prob.bara.subi)); 

% Mosek parameters
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 0.1;
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 0.1;
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 0.1;

%% Start memory logging
% NOTE: This works only on Unix systems. 

hmL = memLog([mfilename,'_',maxcut_data]);
hmL.start;
MEMBEGIN = hmL.prompt;

%% Solve with Mosek

timer = tic;
cputimeBegin = cputime;

[~,res] = mosekopt('minimize info',prob,param);

totalTime = toc(timer);
cputimeEnd = cputime;

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

%% Stop memory logging

MEMEND = hmL.prompt;
hmL.stop;
out.memory = (MEMEND - MEMBEGIN)/1000; %in MB

%% Evaluate errors

res.sol.itr.bars = 'deleted';
res.sol.itr.skc = 'deleted';
res.sol.itr.skx = 'deleted';
res.sol.itr.suc = 'deleted';
res.sol.itr.sux = 'deleted';
res.sol.itr.ux = 'deleted';
res.sol.itr.xx = 'deleted';
res.sol.itr.y = 'deleted';

X = zeros(n);
[ii,jj] = ndgrid(1:n); % ii and jj are row and column indices respectively
X(ii>=jj) = res.sol.itr.barx; % fill in values in column-major 
X = X + tril(X,-1).'; % transpose to get result
res.sol.itr.barx = 'deleted'; % I deleted these for reducing the file size

R = 10; % rank/sketch size parameter
[u,~] = eigs(X, R, 'LM');
cutvalue = 0; 
for t = 1:R
    sign_evec = sign(u(:,t));
    rankvalue = -sign_evec'*(C*sign_evec);
    cutvalue = max(cutvalue, rankvalue);
end

out.cutvalue = cutvalue;
out.primalFeas = norm(diag(X)-b)/(1+norm(b));
out.primalObj = C(:)'*X(:);

%% Save Results

if ~exist(['results/MaxCut/',maxcut_data],'dir'), mkdir(['results/MaxCut/',maxcut_data]); end
save(['results/MaxCut/',maxcut_data,'/Mosek.mat'],'res','out','-v7.3');

%% Last edit: Alp Yurtsever - July 24, 2020
%% Test setup for Primal Dual Convergence (MaxCut SDP) - Solved with Mosek
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
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-3;
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-3;
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-3;

%% Solve with Mosek

timer = tic;
cputimeBegin = cputime;

[~,res] = mosekopt('minimize info',prob,param);

totalTime = toc(timer);
cputimeEnd = cputime;

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

%% Evaluate errors

res.sol.itr.skc = 'deleted';
res.sol.itr.skx = 'deleted';
res.sol.itr.suc = 'deleted';
res.sol.itr.sux = 'deleted';
res.sol.itr.ux = 'deleted';
res.sol.itr.xx = 'deleted';

X = zeros(n);
[ii,jj] = ndgrid(1:n); % ii and jj are row and column indices respectively
X(ii>=jj) = res.sol.itr.barx; % fill in values in column-major 
X = X + tril(X,-1).'; % transpose to get result
res.sol.itr.barx = 'deleted'; % I deleted these for reducing the file size
[ii,jj] = ndgrid(1:n); % ii and jj are row and column indices respectively
Z = zeros(n);
Z(ii>=jj) = res.sol.itr.bars; % fill in values in column-major 
Z = Z + tril(Z,-1).'; % transpose to get result
res.sol.itr.bars = 'deleted';
y = res.sol.itr.y;
res.sol.itr.y = 'deleted';

dobj = -b'*y;
pobj = C(:)'*X(:);
out.dimacs.err1 = norm(diag(X)-b)/(1+norm(b));
out.dimacs.err2 = max(-min(eig(X)),0)/(1+norm(b));
out.dimacs.err3 = norm(diag(y)+Z-C,'fro')/(1+norm(C,'fro'));
out.dimacs.err4 = max(-min(eig(Z)),0)/(1+norm(C,'fro'));
out.dimacs.err5 = (pobj+dobj)/(1+abs(pobj)+abs(dobj));
out.dimacs.err6 = (X(:)'*Z(:))/(1+abs(pobj)+abs(dobj));

R = 10;
[u,~] = eigs(X, R, 'LM');
cutvalue = 0; 
for t = 1:R
    sign_evec = sign(u(:,t));
    rankvalue = -(sign_evec'*(C*sign_evec));
    cutvalue = max(cutvalue, rankvalue);
end

out.cutvalue = cutvalue;
out.primalFeas = norm(diag(X)-b)/(1+norm(b));
out.primalObj = C(:)'*X(:);

%% Save Results

if ~exist(['results/DualMaxCut/',maxcut_data],'dir'), mkdir(['results/DualMaxCut/',maxcut_data]); end
save(['results/DualMaxCut/',maxcut_data,'/Mosek.mat'],'res','out','-v7.3');

%% Last edit: Alp Yurtsever - July 24, 2020
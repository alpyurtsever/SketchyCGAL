%% Test setup for Primal Dual Convergence (MaxCut SDP) - Solved with SeDuMi
%% Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)

%% Choose data
% NOTE: You need to download data from GSET and locate them to under the
% "FilesMaxCut/data/G/" folder (resp. DIMACS10, "FilesMaxCut/data/DIMACS10/").

% maxcut_data = 'G/G1'; % you can choose other data files as well
% maxcut_data = 'DIMACS10/belgium_osm';

%% Preamble
rng(0,'twister');

% Add SEDUMI to path!
% NOTE: You need to download SEDUMI if you don't already have it.

% addpath(genpath('DESTINATION OF SEDUMI'));

% Also make a copy of the 'sedumi.m' as './FilesMaxcut/sedumi_modified/'
% and comment out lines  638-642. See Section SM5.1.4 in our paper for more
% details.
addpath('./FilesMaxcut/sedumi_modified/');

%% load data

load(['./FilesMaxCut/data/',maxcut_data]);

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
C = (-0.25).*L(:);
clearvars L;
indI = (0:n-1)'*(n+1) + 1;
indJ = (1:n)';
valA = ones(n,1);
b = ones(n,1);
At = sparse(indI,indJ,valA,n^2,n);
clearvars indI indJ valA;

pars.eps = 1e-3;

%% Solve with SEDUMI

timer = tic;
cputimeBegin = cputime;

[X,y,info] = sedumi_modified(At,b,C,K,pars);

cputimeEnd = cputime;
totalTime = toc(timer);

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

%% Evaluate errors

C = reshape(C,[n,n]);
X = reshape(X,[n,n]);
X = 0.5*(X+X'); % sedumi does not return symmetric outputs
Z = C - diag(y);

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
save(['results/DualMaxCut/',maxcut_data,'/SeDuMi.mat'],'info','out','-v7.3');

%% Last edit: Alp Yurtsever - July 24, 2020
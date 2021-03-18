%% Test setup for Primal Dual Convergence (MaxCut SDP) - Solved with SDPT3
%% Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)

%% Choose data
% NOTE: You need to download data from GSET and locate them to under the
% "FilesMaxCut/data/G/" folder (resp. DIMACS10, "FilesMaxCut/data/DIMACS10/").

% maxcut_data = 'G/G1'; % you can choose other data files as well
% maxcut_data = 'DIMACS10/belgium_osm';

%% Preamble
rng(0,'twister');

% Add SDPT3 to path!
% NOTE: You need to download SDPT3 if you don't already have it.

% addpath(genpath('DESTINATION OF SDPT3'));

%% load data

load(['./FilesMaxCut/data/',maxcut_data]);

n = size(Problem.A,1);
C{1} = spdiags(Problem.A*ones(n,1),0,n,n) - Problem.A;
clearvars Problem;
C{1} = 0.5*(C{1}+C{1}');
C{1} = -0.25*C{1};
b = ones(n,1);

%% Construct SDPT3 Problem

blk{1,1} = 's';
blk{1,2} = n;
A = cell(1,n);
for k = 1:n; A{k} = sparse(k,k,1,n,n); end
Avec = svec(blk,A,ones(size(blk,1),1));
clearvars A;

OPTIONS.gaptol = 1e-3;

%% Solve with SDPT3

timer = tic;
cputimeBegin = cputime;

[obj,X,y,Z,info,runhist] = sdpt3(blk,Avec,C,b,OPTIONS);

cputimeEnd = cputime;
totalTime = toc(timer);

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

%% Evaluate errors

X = cell2mat(X);
Z = cell2mat(Z);
C = cell2mat(C);

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
save(['results/DualMaxCut/',maxcut_data,'/SDPT3.mat'],'obj','info','runhist','out','-v7.3');

%% Last edit: Alp Yurtsever - July 24, 2020
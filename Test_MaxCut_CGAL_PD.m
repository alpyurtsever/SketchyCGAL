%% Test setup for Primal Dual Convergence (MaxCut SDP) - Solved with SketchyCGAL
%% Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)

%% Choose data
% NOTE: You need to download data from GSET and locate them to under the
% "FilesMaxCut/data/G/" folder (resp. DIMACS10, "FilesMaxCut/data/DIMACS10/").

% maxcut_data = 'G/G1'; % you can choose other data files as well
% maxcut_data = 'DIMACS10/belgium_osm';

%% Preamble
rng(0,'twister');
addpath utils;
addpath solver;

%% Load data

load(['./FilesMaxCut/data/',maxcut_data]);

n = size(Problem.A,1);
C = spdiags(Problem.A*ones(n,1),0,n,n) - Problem.A;
C = 0.5*(C+C'); % symmetrize if not symmetric
C = (-0.25).*C;

clearvars Problem;

%% Construct the Black Box Oracles for MaxCut SDP

Primitive1 = @(x) C*x;
Primitive2 = @(y,x) y.*x;
Primitive3 = @(x) sum(x.^2,2);
a = [0,1.05*n];
b = ones(n,1);

% Compute scaling factors
SCALE_X = 1/n;
SCALE_C = 1/norm(C,'fro');

%% Solve using SketchyCGAL

R = 10; % rank/sketch size parameter
beta0 = 1; % we didn't tune - choose 1 - you can tune this!
K = inf;
maxit = 1e6; % limit on number of iterations

timer = tic;
cputimeBegin = cputime;

[out, U, D, y, AX, pobj] = CGAL( n, Primitive1, Primitive2, Primitive3, a, b, R, maxit, beta0, K, ...
    'FLAG_MULTRANK_P1',true,... % This flag informs that Primitive1 can be applied to find AUU' for any size U. 
    'FLAG_MULTRANK_P3',true,... % This flag informs that Primitive1 can be applied to find (A'y)U for any size U.
    'SCALE_X',SCALE_X,... % SCALE_X prescales the primal variable X of the problem
    'SCALE_C',SCALE_C,... % SCALE_C prescales the cost matrix C of the problem
    'stoptol',1e-3); % algorithm stops when 1e-3 accuracy is achieved
                     
cputimeEnd = cputime;
totalTime = toc(timer);

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

%% Evaluate errors

dobj = b'*y;
Z = C+diag(y);
out.dimacs.err1 = norm(AX-b)/(1+norm(b));
out.dimacs.err2 = 0; % this is theoretically 0 for CGAL
out.dimacs.err3 = 0; % this is theoretically 0 for CGAL
out.dimacs.err4 = max(-min(eig(Z)),0)/(1+norm(C,'fro'));
out.dimacs.err5 = (pobj+dobj)/(1+abs(pobj)+abs(dobj));
out.dimacs.err6 = (pobj + y'*AX)/(1+abs(pobj)+abs(dobj));

cutvalue = 0;
for t = 1:R
    sign_evec = sign(U(:,t));
    rankvalue = -(sign_evec'*(C*sign_evec));
    cutvalue = max(cutvalue, rankvalue);
end

out.cutvalue = cutvalue;
out.primalObj = pobj;
out.primalFeas = norm(AX-b)/(1+norm(b));

%% Save results

if ~exist(['results/DualMaxCut/',maxcut_data],'dir'), mkdir(['results/DualMaxCut/',maxcut_data]); end
save(['results/DualMaxCut/',maxcut_data,'/SketchyCGAL.mat'],'out','-v7.3');

%% Last edit: Alp Yurtsever - July 24, 2020
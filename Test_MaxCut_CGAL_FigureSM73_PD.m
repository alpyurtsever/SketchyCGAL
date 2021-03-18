%% Test setup for Primal Dual Convergence (MaxCut SDP) - Solved with SketchyCGAL (without stopping)
%% Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)

%% Choose data
% NOTE: You need to download data from GSET and locate them to under the
% "FilesMaxCut/data/G/" folder (resp. DIMACS10, "FilesMaxCut/data/DIMACS10/").

maxcut_data = 'G/G67'; % you can choose other data files as well

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

%% Evalaute DIMACS errors
% We skip err2 and err3 because both are provably zero for CGAL

err{1} = 'err1'; err{2} = @(U,AX,pobj,y) err1(C,b,AX,pobj,y);
err{3} = 'err4'; err{4} = @(U,AX,pobj,y) err4();
err{5} = 'err5'; err{6} = @(U,AX,pobj,y) err5();
err{7} = 'err6'; err{8} = @(U,AX,pobj,y) err6();

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
    'EVALSURROGATEGAP',true,... % Evaluates the surrogate gap although we do not use it for stopping
    'errfncs',err); % err defines the spectral rounding for maxcut
                     
cputimeEnd = cputime;
totalTime = toc(timer);

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

out.primalObj = pobj;
out.primalFeas = norm(AX-b)/(1+norm(b));

%% Save results

if ~exist(['results/DualMaxCut/',maxcut_data],'dir'), mkdir(['results/DualMaxCut/',maxcut_data]); end
save(['results/DualMaxCut/',maxcut_data,'/FigureSM73.mat'],'out','-v7.3');

%% Implement DIMACS errors

function out = err1(C,b,AX,pobj,y)
dobj = b'*y;
Z = C+diag(y);
global DIMACS;
DIMACS = {};
DIMACS.err1 = norm(AX-b)/(1+norm(b));
DIMACS.err4 = max(-min(eig(Z)),0)/(1+norm(C,'fro'));
DIMACS.err5 = (pobj+dobj)/(1+abs(pobj)+abs(dobj));
DIMACS.err6 = (pobj + y'*AX)/(1+abs(pobj)+abs(dobj));
out = DIMACS.err1;
end
function out = err4(), global DIMACS; out = DIMACS.err4; end
function out = err5(), global DIMACS; out = DIMACS.err5; end
function out = err6(), global DIMACS; out = DIMACS.err6; end

%% Last edit: Alp Yurtsever - July 24, 2020
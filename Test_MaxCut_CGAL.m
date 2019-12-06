%%  Test setup for MaxCut SDP - Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)
% Solved with SketchyCGAL
% To replicate the results in Figure 7.1 in the paper

%% Choose data
% NOTE: You need to download data from GSET and locate them to under the
% "FilesMaxCut/data/G/" folder (resp. DIMACS10, "FilesMaxCut/data/DIMACS10/").

maxcut_data = 'G/G1'; % you can choose other data files as well
% maxcut_data = 'DIMACS10/belgium_osm';

%% Preamble
rng(0,'twister');
addpath utils;
addpath solver;

%% Load data

% load(['/Users/alp/Documents/MATLAB/DATASETS/MaxCut/',maxcut_data]);
% load(['/home/yurtseve/MATLAB_COMMON/DATASETS/MaxCut/',maxcut_data]);
% load(['./FilesMaxCut/data/',maxcut_data]);

n = size(Problem.A,1);
C = spdiags(Problem.A*ones(n,1),0,n,n) - Problem.A;
C = 0.5*(C+C'); % symmetrize if not symmetric
C = (-0.25).*C;

clearvars Problem;

%% Construct the Black Box Oracles for MaxCut SDP

Primitive1 = @(x) C*x;
Primitive2 = @(y,x) y.*x;
Primitive3 = @(x) sum(x.^2,2);
a = n;
b = ones(n,1);

% Compute scaling factors
SCALE_X = 1/a;
SCALE_C = 1/norm(C,'fro');

%% Implement rounding

err{1} = 'cutvalue'; % name for error
err{2} = @(u) round(C,u); % function definition at the bottom of this script

%% Start memory logging
% NOTE: This works only on Unix systems. 

hmL = memLog([maxcut_data,'cgal']);
hmL.start;
MEMBEGIN = hmL.prompt;

%% Solve using SketchyCGAL

R = 10; % rank/sketch size parameter
beta0 = 1; % we didn't tune - choose 1 - you can tune this!
K = inf;
maxit = 1e6; % limit on number of iterations

timer = tic;
cputimeBegin = cputime;

out = CGAL( n, Primitive1, Primitive2, Primitive3, a, b, R, maxit, beta0, K, ...
    'FLAG_MULTRANK_P1',true,... % err defines the spectral rounding for maxcut
    'FLAG_MULTRANK_P3',true,... % err defines the spectral rounding for maxcut
    'SCALE_X',SCALE_X,... % err defines the spectral rounding for maxcut
    'SCALE_C',SCALE_C,... % err defines the spectral rounding for maxcut
    'errfncs',err,... % err defines the spectral rounding for maxcut
    'stoptol',1e-4); % algorithm stops when 1e-4 accuracy is achieved
                     % NOTE1: stoptol checks the accuracy of the SCALED PROBLEM!
                     % NOTE2: remove 'stoptol' option to run method for 1e6 itearions to generate the results in Figure 7.2

cputimeEnd = cputime;
totalTime = toc(timer);

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

out.cutvalue_r = out.info.cutvalue(end); % out.cutvalue_r = round_r(C,U);
out.primalObj = out.info.primalObj(end); % out.primalObj = trace(U'*C*U*D);
out.primalFeas = out.info.primalFeas(end); % out.primalFeas = norm(diag(U*D*U') - b)/max(norm(b),1);
out.sdpFeas = 0; % abs(min(min(diag(D)),0)); % should be always 0 for CGAL

%% Stop memory logging

MEMEND = hmL.prompt;
hmL.stop;
out.memory = (MEMEND - MEMBEGIN)/1000; %in MB

%% Save results

if ~exist(['results/MaxCut/',maxcut_data],'dir'), mkdir(['results/MaxCut/',maxcut_data]); end
save(['results/MaxCut/',maxcut_data,'/SketchyCGAL.mat'],'out','-v7.3');

%% Implement rounding
% NOTE: Defining a function in a MATLAB script was not available in older
% versions. If you are using an old version of MATLAB, you might need to
% save this function as a seperate ".m" file in your path.

function cutvalue = round(C,u)
cutvalue = 0;
for t = 1:size(u,2)
    sign_evec = sign(u(:,t));
    rankvalue = -(sign_evec'*(C*sign_evec));
    cutvalue = max(cutvalue, rankvalue);
end
end

%% Last edit: Alp Yurtsever - November 18, 2019
%%  Test setup for MaxCut SDP - Solved with SketchyCGAL (without stopping)
%% Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)

%% Choose data
% NOTE: You need to download data from GSET and locate them to under the
% "FilesMaxCut/data/G/" folder (resp. DIMACS10, "FilesMaxCut/data/DIMACS10/").

R = 10; % rank/sketch size parameter {10,25,100}
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
a = n;
b = ones(n,1);

SCALE_X = 1/n;
SCALE_C = 1/norm(C,'fro');

%% Implement rounding

err{1} = 'cutvalue'; % name for error
err{2} = @(u) round(C,u); % function definition at the bottom of this script

%% Solve using SketchyCGAL

beta0 = 1; % we didn't tune - choose 1 - you can tune this!
K = inf;
maxit = 1e6; % limit on number of iterations

out = CGAL( n, Primitive1, Primitive2, Primitive3, a, b, R, maxit, beta0, K, ...
    'FLAG_MULTRANK_P1',true,... % This flag informs that Primitive1 can be applied to find AUU' for any size U. 
    'FLAG_MULTRANK_P3',true,... % This flag informs that Primitive1 can be applied to find (A'y)U for any size U.
    'SCALE_X',SCALE_X,... % SCALE_X prescales the primal variable X of the problem
    'SCALE_C',SCALE_C,... % SCALE_C prescales the cost matrix C of the problem
    'errfncs',err); % err defines the spectral rounding for maxcut

%% Save results

if ~exist(['results/MaxCut/',maxcut_data],'dir'), mkdir(['results/MaxCut/',maxcut_data]); end
save(['results/MaxCut/',maxcut_data,'/Figure71-R=',num2str(R),'.mat'],'out','-v7.3');

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

%% Last edit: Alp Yurtsever - July 24, 2020
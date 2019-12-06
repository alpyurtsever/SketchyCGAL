%% Choose data
% NOTE: You need to download data from GSET and locate them to under the
% "FilesMaxCut/data/G/" folder (resp. DIMACS10, "FilesMaxCut/data/DIMACS10/").

dataNo = 1; % choose from 1 to 10

%% Preamble
rng(0,'twister');
addpath solver;

% Add manopt to path!

% addpath(genpath('DESTINATION OF MANOPT'));
addpath(genpath('./FilesMaxcut/manopt'));

%% Load data

load(['./FilesMaxCut/data/WALDSPURGER/C',num2str(dataNo),'.mat']);
n = size(C,1);
optval = xopt'*C*xopt;

%% Search rank - sketch size

% Note that the solution has rank 1
R = 2;

%% Solve with Manopt

% The fixed rank elliptope geometry describes symmetric, positive
% semidefinite matrices of size n with rank r and all diagonal entries
% are 1.

% Set options
options.verbosity = 1;
options.maxiter = 1e3;
options.tolgradnorm = -inf;

% with rank 2
manifold = elliptopefactory(n, 2);
problem.M = manifold;
problem.cost  = @(Y) trace(Y'*C*Y);
problem.egrad = @(Y) 2*C*Y;
problem.ehess = @(Y, U) 2*(C*U);

[~, ~, infoRank2] = trustregions(problem,[],options);

%% Solve with SketchyCGAL

Primitive1 = @(x) C*x;
Primitive2 = @(y,x) y.*x;
Primitive3 = @(x) sum(x.^2,2);
a = n;
b = ones(n,1);

beta0 = 1; % we didn't tune - choose 1 - you can tune this!
K = inf;
maxit = 1e5; % limit on number of iterations

out = CGAL( n, Primitive1, Primitive2, Primitive3, a, b, R, maxit, beta0, K);

%% Plots
close all

hfig1 = figure('Position',[100,100,1200,260]);
set(hfig1,'name','MaxCut-Waldspurger','numbertitle','off');

subplot(141)
loglog([infoRank2.iter], [infoRank2.gradnorm],'Color',[0.75,0,0]);
ylabel('gradient norm','Interpreter','latex','FontSize',17);
ht = title('manopt');
ht.Interpreter = 'latex';
ylim([1e-13,1e4])
xlim([1,1e3])
set(gca,'YTick',10.^(-21:3:20));

subplot(142)
loglog([infoRank2.iter], abs([infoRank2.cost]-optval)/max(1,abs(optval)),'Color',[0.75,0,0]);
ylim([1e1,1e4])
xlim([1,1e3])
ylabel('objective residual','Interpreter','latex','FontSize',17);
ht = title('manopt');
ht.Interpreter = 'latex';
hl = legend('$R = 2$');
hl.Interpreter = 'latex';
hl.Location = 'SouthWest';
hl.FontSize = 14;
set(gca,'YTick',10.^(-10:10));

subplot(143)
loglog(out.iteration,abs(out.info.primalObj-optval)/max(1,abs(optval)),'Color',[0,0,0.75],'LineStyle','--');
hold on
ylim([0.999e-5,1e2])
set(gca,'YTick',10.^(-10:10));
loglog(out.iteration,abs(out.info.skPrimalObj-optval)/max(1,abs(optval)),'Color',[0,0,0.75]);
ylabel('objective residual','Interpreter','latex','FontSize',17);
ht = title('SketchyCGAL');
ht.Interpreter = 'latex';

subplot(144)
hl1 = loglog(out.iteration,out.info.primalFeas,'Color',[0,0,0.75],'LineStyle','--');
hold on
hl2 = loglog(out.iteration,out.info.skPrimalFeas,'Color',[0,0,0.75]);
ylim([1e-10,1])
ylabel('infeasibility','Interpreter','latex','FontSize',17);
ht = title('SketchyCGAL');
ht.Interpreter = 'latex';
set(gca,'YTick',10.^(-10:2:10));

hl = legend([hl2,hl1], '$R = 2$','Implicit');
hl.Interpreter = 'latex';
hl.Location = 'SouthWest';
hl.FontSize = 14;

for t = 1:4
    subplot(1,4,t)
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
    set(gca,'FontSize',14,'TickLabelInterpreter','latex');
    grid on, grid minor, grid minor;
    set(gca,'XTick',10.^(0:10));
    xlabel('iteration','Interpreter','latex','FontSize',16);
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
end

drawnow;

return;
%% Pass to second experiment
%% Get statistics about how many times manopt converges to spurious critical point
clearvars;
rng(0,'twister');
% addpath(genpath('DESTINATION OF MANOPT'));
addpath(genpath('./FilesMaxcut/manopt'));

n = 100;
Rmax = floor(sqrt(2*n+9/4)-3/2)+1; % factorization rank
stats = nan(10,length(Rmax));

options.verbosity = 0;

for dataNo = 1:10
    load(['./FilesMaxcut/data/WALDSPURGER/C',num2str(dataNo),'.mat']);
    for R = 2:Rmax
        numFail = 0;
        for t = 1:100
            
            manifold = elliptopefactory(n, R);
            problem.M = manifold;
            problem.cost  = @(Y) trace(Y'*C*Y);
            problem.egrad = @(Y) 2*C*Y;
            problem.ehess = @(Y, U) 2*(C*U);
            
            [~, Ycost] = trustregions(problem,[],options);
            if Ycost > 1e-3
                numFail = numFail + 1;
            end
            
        end
        stats(dataNo,R-1) = numFail;
        fprintf('Data C%d: manopt with R = %d failed at %d or 100 instances.\n',dataNo,R,numFail);
    end
end

%% Save results

if ~exist('results/MaxCut/Waldspurger/','dir'), mkdir('results/MaxCut/Waldspurger/'); end
save('results/MaxCut/Waldspurger/stats.mat','stats');

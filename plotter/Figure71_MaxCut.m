%% Script used to generate Figure 7.1
% Before running this, you should run the MaxCut experiment with G67 
% dataset with R = {10,25,100}. This script only generates the figure from 
% the saved results of these experiments. 

%% Preamble
close all; clearvars;

%% Load results
data.CGAL10 = load('../results/MaxCut/G/G67/Figure71-R=10.mat');
data.CGAL25 = load('../results/MaxCut/G/G67/Figure71-R=25.mat');
data.CGAL100 = load('../results/MaxCut/G/G67/Figure71-R=100.mat');
data.SDPT3 = load('../results/MaxCut/G/G67/SDPT3-GroundTruth.mat');
optval = data.SDPT3.obj(1); % Ground truth approximated with SDPT3

%% Plot Fig1: infeasibility vs iteration counter
hfig1 = figure('Position',[100,100,400,260]);
set(hfig1 ,'name','G67-infeasibility-vs-iteration','numbertitle','off');

loglog(data.CGAL10.out.iteration,   data.CGAL10.out.info.skPrimalFeas); hold on;
loglog(data.CGAL25.out.iteration,   data.CGAL25.out.info.skPrimalFeas);
loglog(data.CGAL100.out.iteration,  data.CGAL100.out.info.skPrimalFeas);
loglog(data.CGAL10.out.iteration,   data.CGAL10.out.info.primalFeas,'k');

ylabel('infeasibility','Interpreter','latex','FontSize',16);
xlabel('iteration','Interpreter','latex','FontSize',16);

axis tight
set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
ax = gca;
ax.XTick = 10.^(-10:10);
ax.YTick = 10.^(-10:10);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
ax.Box = 'on';
set(findall(ax, 'Type', 'Line'),'LineWidth',2);
grid on, grid minor, grid minor;

%% Plot Fig2: objective residual vs iteration counter
hfig2 = figure('Position',[100,100,400,260]);
set(hfig2 ,'name','G67-objective-vs-iteration','numbertitle','off');

loglog(data.CGAL10.out.iteration,   abs(data.CGAL10.out.info.skPrimalObj - optval)/max(1,abs(optval))); hold on;
loglog(data.CGAL25.out.iteration,   abs(data.CGAL25.out.info.skPrimalObj - optval)/max(1,abs(optval)));
loglog(data.CGAL100.out.iteration,  abs(data.CGAL100.out.info.skPrimalObj - optval)/max(1,abs(optval)));
loglog(data.CGAL10.out.iteration,   abs(data.CGAL10.out.info.primalObj - optval)/max(1,abs(optval)),'k');

ylabel('objective residual','Interpreter','latex','FontSize',16);
xlabel('iteration','Interpreter','latex','FontSize',16);

axis tight
set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
ax = gca;
ax.XTick = 10.^(-10:10);
ax.YTick = 10.^(-10:10);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
ax.Box = 'on';
set(findall(ax, 'Type', 'Line'),'LineWidth',2);
grid on, grid minor, grid minor;

hl = legend('$R = 10$','$R = 25$','$R = 100$','Full matrix (implicit)');
hl.Interpreter = 'latex';
hl.FontSize = 13;
hl.Location = 'SouthWest';

%% Plot Fig3: weight of cut vs iteration counter
hfig3 = figure('Position',[100,100,400,260]);
set(hfig3 ,'name','G67-weight-vs-iteration','numbertitle','off');

semilogx(data.CGAL10.out.iteration,     data.CGAL10.out.info.cutvalue); hold on;
semilogx(data.CGAL25.out.iteration,     data.CGAL25.out.info.cutvalue);
semilogx(data.CGAL100.out.iteration,    data.CGAL100.out.info.cutvalue);
semilogx(data.CGAL10.out.iteration,     data.SDPT3.out.cutvalue*ones(size(data.CGAL10.out.iteration)),'LineStyle','--','Color','black');

ylabel('weight of cut','Interpreter','latex','FontSize',16);
xlabel('iteration','Interpreter','latex','FontSize',16);

axis tight
yLim = ylim; yLim(2) = 6500; ylim(yLim);
set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
ax = gca;
ax.XTick = 10.^(-10:10);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
ax.Box = 'on';
ax.YRuler.MinorTick = 'off'; %or 'off'
set(findall(ax, 'Type', 'Line'),'LineWidth',2);
grid on, grid minor, grid minor;

%% Plot Fig4: infeasibility vs time
hfig4 = figure('Position',[100,100,400,260]);
set(hfig4 ,'name','G67-infeasibility-vs-time','numbertitle','off');

loglog(data.CGAL10.out.time,     data.CGAL10.out.info.skPrimalFeas); hold on;
loglog(data.CGAL25.out.time,     data.CGAL25.out.info.skPrimalFeas);
loglog(data.CGAL100.out.time,    data.CGAL100.out.info.skPrimalFeas);
loglog(data.CGAL10.out.time,     data.CGAL10.out.info.primalFeas,'k');

ylabel('infeasibility','Interpreter','latex','FontSize',16);
xlabel('time (sec)','Interpreter','latex','FontSize',16);

axis tight
set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
ax = gca;
ax.XTick = 10.^(-10:10);
ax.YTick = 10.^(-10:10);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
ax.Box = 'on';
set(findall(ax, 'Type', 'Line'),'LineWidth',2);
grid on, grid minor, grid minor;

%% Plot Fig5: objective suboptimality vs time
hfig5 = figure('Position',[100,100,400,260]);
set(hfig5 ,'name','G67-objective-vs-time','numbertitle','off');

loglog(data.CGAL10.out.time,     abs(data.CGAL10.out.info.skPrimalObj - optval)/max(1,abs(optval))); hold on;
loglog(data.CGAL25.out.time,     abs(data.CGAL25.out.info.skPrimalObj - optval)/max(1,abs(optval)));
loglog(data.CGAL100.out.time,    abs(data.CGAL100.out.info.skPrimalObj - optval)/max(1,abs(optval)));
loglog(data.CGAL10.out.time,     abs(data.CGAL10.out.info.primalObj - optval)/max(1,abs(optval)),'k');

ylabel('objective residual','Interpreter','latex','FontSize',16);
xlabel('time (sec)','Interpreter','latex','FontSize',16);

axis tight
set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
ax = gca;
ax.XTick = 10.^(-10:10);
ax.YTick = 10.^(-10:10);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
ax.Box = 'on';
set(findall(ax, 'Type', 'Line'),'LineWidth',2);
grid on, grid minor, grid minor;

%% Plot Fig6: weight of cut vs time
hfig6 = figure('Position',[100,100,400,260]);
set(hfig6 ,'name','G67-weight-vs-time','numbertitle','off');

semilogx(data.CGAL10.out.time,   data.CGAL10.out.info.cutvalue); hold on;
semilogx(data.CGAL25.out.time,   data.CGAL25.out.info.cutvalue);
semilogx(data.CGAL100.out.time,  data.CGAL100.out.info.cutvalue);
semilogx(data.CGAL10.out.time,   data.SDPT3.out.cutvalue*ones(size(data.CGAL10.out.time)),'k--');

ylabel('weight of cut','Interpreter','latex','FontSize',16);
xlabel('time (sec)','Interpreter','latex','FontSize',16);

axis tight
yLim = ylim; yLim(2) = 6500; ylim(yLim);
set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
ax = gca;
ax.XTick = 10.^(-10:10);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
ax.Box = 'on';
ax.YRuler.MinorTick = 'off'; %or 'off'
set(findall(ax, 'Type', 'Line'),'LineWidth',2);
grid on, grid minor, grid minor;

%% Last edit: Alp Yurtsever - July 24, 2020

%% Script used to generate Figure SM7.3
% Before running this, you should run the MaxCut experiment with G67 
% dataset with R = {10,25,100}. This script only generates the figure from 
% the saved results of these experiments. 

%% Preamble
close all; clearvars;

%% Load results
data = load('../results/DualMaxCut/G/G67/FigureSM73.mat');

%% Plot Fig1: infeasibility vs iteration counter
hfig1 = figure('Position',[100,100,400,260]);
set(hfig1 ,'name','G67-dimacs-err1','numbertitle','off');

loglog(data.out.iteration, data.out.info.err1,'k'); hold on;

ylabel('err$_1$','Interpreter','latex','FontSize',16);
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
set(hfig2 ,'name','G67-dimacs-err4','numbertitle','off');

loglog(data.out.iteration, data.out.info.err4,'k'); hold on;

ylabel('err$_4$','Interpreter','latex','FontSize',16);
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

%% Plot Fig3: weight of cut vs iteration counter
hfig3 = figure('Position',[100,100,400,260]);
set(hfig3 ,'name','G67-dimacs-err5','numbertitle','off');

loglog(data.out.iteration, abs(data.out.info.err5),'k'); hold on;

ylabel('$|$err$_5|$','Interpreter','latex','FontSize',16);
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

%% Plot Fig4: infeasibility vs cpu time
hfig4 = figure('Position',[100,100,400,260]);
set(hfig4 ,'name','G67-dimacs-err6','numbertitle','off');

loglog(data.out.iteration, abs(data.out.info.err6),'k'); hold on;

ylabel('$|$err$_6|$','Interpreter','latex','FontSize',16);
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

%% Plot Fig3: weight of cut vs iteration counter
hfig5 = figure('Position',[100,100,400,260]);
set(hfig5 ,'name','G67-surrogate-gap','numbertitle','off');

loglog(data.out.iteration, abs(data.out.info.stopObj),'k'); hold on;

ylabel('$|$surrogate~gap$|$','Interpreter','latex','FontSize',16);
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

%% Last edit: Alp Yurtsever - July 24, 2020

%% Script used to generate Figure 7.3 
% Before running this, you should run the PhaseRetrieval experiment with  
% all methods for all problem sizes (and Monte Carlo iterations). This 
% script only generates the figure from the saved results of these 
% experiments.

clearvars;
close all;
nSweep = [1e2,1e3,1e4,1e5,1e6];

warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');

CGAL.mem = nan(length(nSweep),20);
CGAL.time = nan(length(nSweep),20);
CGAL.cputime = nan(length(nSweep),20);
ptr = 0;
for n = nSweep
    ptr = ptr + 1;
    for MC = 1:20
        if exist(['../results/PhaseRetrieval/10/',num2str(n),'_',num2str(MC),'_CGAL.mat'],'file')
            data = load(['../results/PhaseRetrieval/10/',num2str(n),'_',num2str(MC),'_CGAL.mat']);
            CGAL.time(ptr,MC) = data.out.totalTime;
            CGAL.mem(ptr,MC) = data.out.memory;
            CGAL.cputime(ptr,MC) = data.out.totalCpuTime;
            clearvars data;
        else
            CGAL.mem(ptr,MC) = nan;
            CGAL.time(ptr,MC) = nan;
            CGAL.cputime(ptr,MC) = nan;
        end
    end
end

ThinCGAL.mem = nan(length(nSweep),20);
ThinCGAL.time = nan(length(nSweep),20);
ThinCGAL.cputime = nan(length(nSweep),20);
ptr = 0;
for n = nSweep
    ptr = ptr + 1;
    for MC = 1:20
        if exist(['../results/PhaseRetrieval/10/',num2str(n),'_',num2str(MC),'_ThinCGAL.mat'],'file')
            data = load(['../results/PhaseRetrieval/10/',num2str(n),'_',num2str(MC),'_ThinCGAL.mat']);
            ThinCGAL.time(ptr,MC) = data.out.totalTime;
            ThinCGAL.mem(ptr,MC) = data.out.memory;
            ThinCGAL.cputime(ptr,MC) = data.out.totalCpuTime;
            clearvars data;
        else
            ThinCGAL.mem(ptr,MC) = nan;
            ThinCGAL.time(ptr,MC) = nan;
            ThinCGAL.cputime(ptr,MC) = nan;
        end
    end
end

SketchyCGAL.mem = nan(length(nSweep),20);
SketchyCGAL.time = nan(length(nSweep),20);
SketchyCGAL.cputime = nan(length(nSweep),20);
ptr = 0;
for n = nSweep
    ptr = ptr + 1;
    for MC = 1:20
        if exist(['../results/PhaseRetrieval/10/',num2str(n),'_',num2str(MC),'_SketchyCGAL.mat'],'file')
            data = load(['../results/PhaseRetrieval/10/',num2str(n),'_',num2str(MC),'_SketchyCGAL.mat']);
            SketchyCGAL.time(ptr,MC) = data.out.totalTime;
            SketchyCGAL.mem(ptr,MC) = data.out.memory;
            SketchyCGAL.cputime(ptr,MC) = data.out.totalCpuTime;
            clearvars data;
        else
            SketchyCGAL.mem(ptr,MC) = nan;
            SketchyCGAL.time(ptr,MC) = nan;
            SketchyCGAL.cputime(ptr,MC) = nan;
        end
    end
end

SketchyCGAL.mean.mem = mean(SketchyCGAL.mem,2,'omitnan');
CGAL.mean.mem = mean(CGAL.mem,2,'omitnan');
ThinCGAL.mean.mem = mean(ThinCGAL.mem,2,'omitnan');

SketchyCGAL.min.mem = SketchyCGAL.mean.mem - min(SketchyCGAL.mem,[],2,'omitnan');
CGAL.min.mem = CGAL.mean.mem - min(CGAL.mem,[],2,'omitnan');
ThinCGAL.min.mem = ThinCGAL.mean.mem - min(ThinCGAL.mem,[],2,'omitnan');
SketchyCGAL.max.mem = max(SketchyCGAL.mem,[],2,'omitnan') - SketchyCGAL.mean.mem;
CGAL.max.mem = max(CGAL.mem,[],2,'omitnan') - CGAL.mean.mem;
ThinCGAL.max.mem = max(ThinCGAL.mem,[],2,'omitnan') - ThinCGAL.mean.mem;

%% Plot the graph
hfig1 = figure('Position',[100,100,500,300]);
set(hfig1,'name','PhaseRetrieval-StorageVsSize','numbertitle','off');

Y = [CGAL.mean.mem, ThinCGAL.mean.mem, SketchyCGAL.mean.mem];
bar(Y,'grouped');
hold on
errorbar((1:5)-0.225,CGAL.mean.mem,CGAL.min.mem,CGAL.max.mem,'LineWidth',1.5,'LineStyle','none','Color','black');
hold on
errorbar((1:5),ThinCGAL.mean.mem,ThinCGAL.min.mem,ThinCGAL.max.mem,'LineWidth',1.5,'LineStyle','none','Color','black');
errorbar((1:5)+0.225,SketchyCGAL.mean.mem,SketchyCGAL.min.mem,SketchyCGAL.max.mem,'LineWidth',1.5,'LineStyle','none','Color','black');

hl = legend('CGAL','ThinCGAL','SketchyCGAL');
hl.Interpreter = 'latex';
hl.FontSize = 14;
hl.Location = 'NorthWest';

ax = gca;
ax.YScale = 'log';
grid on
axis tight;
ylim([4,2^14])
set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);

ax.XTick = 1:5;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ax.YTick = [4,16,64,256,1024,4096,16384];
ax.XTickLabel = {'$10^2$','$10^3$','$10^4$','$10^5$','$10^6$'};
ax.YTickLabel = {'4 MB','16 MB','64MB','256MB','1GB','4GB','16GB'};
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
ax.Box = 'on';

ax.YRuler.MinorTick = 'off'; %or 'off'
ax.XRuler.MinorTick = 'off'; %or 'off'

ylabel('storage','Interpreter','latex','FontSize',16);
xlabel('problem size: $n$','Interpreter','latex','FontSize',16);

%%
SketchyCGAL.mean.cputime = mean(SketchyCGAL.cputime,2,'omitnan');
CGAL.mean.cputime = mean(CGAL.cputime,2,'omitnan');
ThinCGAL.mean.cputime = mean(ThinCGAL.cputime,2,'omitnan');

SketchyCGAL.min.cputime = SketchyCGAL.mean.cputime - min(SketchyCGAL.cputime,[],2,'omitnan');
CGAL.min.cputime = CGAL.mean.cputime - min(CGAL.cputime,[],2,'omitnan');
ThinCGAL.min.cputime = ThinCGAL.mean.cputime - min(ThinCGAL.cputime,[],2,'omitnan');
SketchyCGAL.max.cputime = max(SketchyCGAL.cputime,[],2,'omitnan') - SketchyCGAL.mean.cputime;
CGAL.max.cputime = max(CGAL.cputime,[],2,'omitnan') - CGAL.mean.cputime;
ThinCGAL.max.cputime = max(ThinCGAL.cputime,[],2,'omitnan') - ThinCGAL.mean.cputime;

hfig2 = figure('Position',[100,100,500,300]);
set(hfig2,'name','PhaseRetrieval-TimeVsSize','numbertitle','off');

YY = [CGAL.mean.cputime, ThinCGAL.mean.cputime, SketchyCGAL.mean.cputime];
bar(YY,'grouped');
hold on
errorbar((1:5)-0.225,CGAL.mean.cputime,CGAL.min.cputime,CGAL.max.cputime,'LineWidth',1.5,'LineStyle','none','Color','black');
hold on
errorbar((1:5),ThinCGAL.mean.cputime,ThinCGAL.min.cputime,ThinCGAL.max.cputime,'LineWidth',1.5,'LineStyle','none','Color','black');
errorbar((1:5)+0.225,SketchyCGAL.mean.cputime,SketchyCGAL.min.cputime,SketchyCGAL.max.cputime,'LineWidth',1.5,'LineStyle','none','Color','black');

ax = gca;
ax.YScale = 'log';
grid on
axis tight;

ax.YTick = 10.^(-100:100);
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ax.XTickLabel = {'$10^2$','$10^3$','$10^4$','$10^5$','$10^6$'};
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
ax.Box = 'on';
ax.YRuler.MinorTick = 'off'; 
ax.XRuler.MinorTick = 'off';

set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);

ylabel('cpu time (sec)','Interpreter','latex','FontSize',16);
xlabel('problem size: $n$','Interpreter','latex','FontSize',16);

%% Last edit: Alp Yurtsever - December 2, 2019

%% Script used to generate Figure 7.5
% Before running this, you should run the QAP experiments with all
% datasets. This script only generates the figure from the saved results 
% of these  experiments.

clearvars
close all
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');

% Table content is from 
% "Semidefinite Programming Approach for the Quadratic Assignment Problem
% with a Sparse Graph", Jose F. S. Bravo Ferreira, Yuehaw Khoo, Amit Singer
% TABLE 4 - https://arxiv.org/pdf/1703.09339.pdf

QAPLIB = {...
    'chr12a',   9552,    34.5,   6.0,    0.0,    42.7    ;...
    'chr12b',   9742,    38.9,   25.4,   11.9,   38.1    ;...
    'chr12c',   11156,   5.8,    2.3,    2.3,    18.6    ;...
    'chr15a',   9896,    2.1,    2.1,    2.1,    52.0    ;...
    'chr15b',   7990,    26.3,   34.5,   29.2,   158.6   ;...
    'chr15c',   9504,    0.0,    0.0,    0.0,    63.3    ;...
    'chr18a',   11098,   69.8,   0.2,    0.2,    76.3    ;...
    'chr18b',   1534,    8.9,    22.9,   29.5,   99.3    ;...
    'chr20a',   2192,    122.5,  76.1,   43.8,   95.4    ;...
    'chr20b',   2298,    62.9,   9.3,    9.3,    82.2    ;...
    'chr20c',   14142,   173.0,  100.1,  111.5,  88.9    ;...
    'chr22a',   6156,    17.2,   7.6,    3.0,    38.3    ;...
    'chr22b',   6194,    7.3,    2.3,    1.0,    40.4    ;...
    'chr25a',   3796,    107.0,  49.2,   25.2,   69.9    ;...
    'esc16a',   68,      8.8,    11.8,   11.8,   11.8    ;...
    'esc16b',   292,     0.0,    0.7,    0.0,    2.7     ;...
    'esc16c',   160,     5.0,    7.5,    8.7,    6.3     ;...
    'esc16d',   16,      12.5,   50.0,   25.0,   75.0    ;...
    'esc16e',   28,      14.3,   7.1,    14.3,   21.4    ;...
    'esc16g',   26,      7.7,    0.0,    15.4,   15.4    ;...
    'esc16h',   996,     1.6,    0.0,    1.6,    16.9    ;...
    'esc16i',   14,      0.0,    0.0,    0.0,    57.1    ;...
    'esc16j',   8,       0.0,    0.0,    0.0,    75.0    ;...
    'esc32a',   130,     115.4,  124.6,  113.8,  93.8    ;...
    'esc32b',   168,     109.5,  114.3,  111.9,  88.1    ;...
    'esc32c',   642,     12.8,   15.9,   13.7,   7.8     ;...
    'esc32d',   200,     38.0,   36.0,   39.0,   21.0    ;...
    'esc32e',   2,       0.0,    0.0,    0.0,    600.0   ;...
    'esc32g',   6,       0.0,    0.0,    0.0,    366.7   ;...
    'esc32h',   438,     24.7,   26.9,   22.8,   18.3    ;...
    'esc64a',   116,     60.3,   53.4,   60.3,   106.9   ;...
    'esc128',   64,      250.0,  206.3,  175.0,  221.9   ;...
    'ste36a',   9526,    70.2,   74.7,   74.2,   76.3     ;...
    'ste36b',   15852,   188.8,  204.3,  211.9,  158.6   ;...
    'ste36c',   8239110, 66.0,   62.8,   63.7,   83.2    ...
    };

% Table content is from 
% "Semidefinite Programming Approach for the Quadratic Assignment Problem
% with a Sparse Graph", Jose F. S. Bravo Ferreira, Yuehaw Khoo, Amit Singer
% TABLE 6 - https://arxiv.org/pdf/1703.09339.pdf

TSPLIB = {...
    'att48',        10628,  213.0,  236.5,  233.6,  329.8   ;...
    'bayg29',       1610,   114.3,  115.8,  114.3,  210.1   ;...
    'bays29',       2020,   107.6,  118.3,  115.4,  164.8   ;...
    'berlin52',     7542,   175.0,  127.2,  127.2,  280.6   ;...
    'bier127',      118282, 216.4,  193.8,  193.8,  234.2   ;...
    'brazil58',     25395,  248.0,  200.8,  200.8,  337.0   ;...
    'burma14',      3323,   24.6,   28.4,   32.3,   95.5    ;...
    'ch130',        6110,   352.4,  380.6,  380.6,  621.3   ;...
    'ch150',        6528,   346.9,  318.2,  318.2,  689.3   ;...
    'dantzig42',    699,    193.1,  174.0,  174.0,  82.0    ;...
    'eil101',       629,    227.3,  235.3,  235.3,  437.7   ;...
    'eil51',        426,    203.6,  205.4,  205.5,  244.4   ;...
    'eil76',        538,    282.9,  183.0,  183.0,  328.2   ;...
    'fri26',        937,    91.6,   39.4,   39.4,   41.6    ;...
    'gr120',        6942,   445.2,  261.6,  261.6,  617.6   ;...
    'gr137',        69853,  264.6,  220.3,  220.3,  38.9    ;...
    'gr17',         2085,   46.8,   32.4,   44.9,   86.9    ;...
    'gr21',         2707,   94.5,   69.7,   66.3,   185.7   ;...
    'gr24',         1272,   89.2,   86.2,   73.9,   129.4   ;...
    'gr48',         5046,   210.2,  187.4,  187.4,  270.4   ;...
    'gr96',         55209,  228.9,  201.7,  201.7,  46.0    ;...
    'hk48',         11461,  222.4,  207.7,  207.7,  281.6   ;...
    'kroA100',      21282,  469.6,  469.0,  469.0,  720.2   ;...
    'kroA150',      26524,  411.0,  467.4,  467.4,  945.8   ;...
    'kroB100',      22141,  411.9,  313.6,  313.6,  624.2   ;...
    'kroB150',      26130,  417.3,  353.7,  353.7,  844.7   ;...
    'kroC100',      20749,  507.4,  445.1,  445.1,  763.0   ;...
    'kroD100',      21294,  504.2,  349.8,  349.8,  654.4   ;...
    'kroE100',      22068,  489.5,  346.3,  346.3,  684.2   ;...
    'lin105',       14379,  303.1,  234.8,  234.8,  248.4   ;...
    'pr107',        44303,  181.5,  207.9,  207.9,  41.6    ;...
    'pr124',        59030,  293.8,  180.2,  180.2   67.6    ;...
    'pr136',        96772,  325.5,  164.7,  164.7,  196.6   ;...
    'pr144',        58537,  255.0,  283.7,  283.7,  59.8    ;...
    'pr76',         108159, 192.2,  194.0,  194.0,  39.4    ;...
    'rat99',        1211,   236.4,  161.5,  161.5,  444.1   ;...
    'rd100',        7910,   438.4,  375.3,  375.3,  506.5   ;...
    'st70',         675,    300.9,  320.0,  317.9,  387.9   ;...
    'swiss42',      1273,   163.2,  190.4,  190.8,  194.0   ;...
    'ulysses16',    6859,   23.6,   20.2,   23.2,   82.7    ;...
    'ulysses22',    7013,   64.5,   57.0,   59.7,   126.3   ...
    };

%% QAP

Optimals = cell2mat(QAPLIB(:,2));

for t = 1:size(QAPLIB,1)
    dataName = QAPLIB{t,1};
    
    load(['../results/QAP/qapdata/',dataName,'/SketchyCGAL.mat']);
    k = sum(~isnan(out.info.eigsRound));
    CGAL.UpperBound(t,1) = min(out.info.eigsRound(1:k));
    CGAL.Iteration(t,1) = out.iteration(k);
    CGAL.Feasibility(t,1) = out.info.primalFeas(k);
    CGAL.Objective(t,1) = out.info.primalObj(k);
end

Gap.CGAL = (CGAL.UpperBound - Optimals)./Optimals * 100;

Gap.CSDP2 = cell2mat(QAPLIB(:,3));
Gap.CSDP3 = cell2mat(QAPLIB(:,4));
Gap.CSDP4 = cell2mat(QAPLIB(:,5));
Gap.PATH = cell2mat(QAPLIB(:,6));

hfig1 = figure('Position',[100,100,1200,300]);
set(hfig1,'name','QAPLIB-barchart','numbertitle','off');

Y = [Gap.CGAL, min([Gap.CSDP2, Gap.CSDP3, Gap.CSDP4],[],2), Gap.PATH];
bar(Y,'grouped');
ax1 = gca;
ax1.XTick = 1:size(QAPLIB,1);
ax1.YTick = 0:100:1000;
set(ax1,'xticklabel',QAPLIB(:,1))
ylabel('relative gap \%','interpreter','latex','fontsize',13)

ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 11;
xtickangle(45)
set(ax1, 'YGrid', 'on', 'XGrid', 'off')
set(gca,'LineWidth',0.5,'TickLength',[0 0]);
hl = legend({'SketchyCGAL','CSDP','PATH'});
hl.Interpreter = 'latex';
hl.FontSize = 12;


%% TSP

Optimals = cell2mat(TSPLIB(:,2));

for t = 1:size(TSPLIB,1)
    dataName = TSPLIB{t,1};
    
    load(['../results/QAP/tspdata/',dataName,'/SketchyCGAL.mat']);
    k = sum(~isnan(out.info.eigsRound));
    CGAL.UpperBound(t,1) = min(out.info.eigsRound(1:k));
    CGAL.Iteration(t,1) = out.iteration(k);
    CGAL.Feasibility(t,1) = out.info.primalFeas(k);
    CGAL.Objective(t,1) = out.info.primalObj(k);
end

Gap.CGAL = (CGAL.UpperBound - Optimals)./Optimals * 100;

Gap.CSDP2 = cell2mat(TSPLIB(:,3));
Gap.CSDP3 = cell2mat(TSPLIB(:,4));
Gap.CSDP4 = cell2mat(TSPLIB(:,5));
Gap.PATH = cell2mat(TSPLIB(:,6));

hfig2 = figure('Position',[100,100,1200,300]);
set(hfig2,'name','TSPLIB-barchart','numbertitle','off');
Y = [Gap.CGAL, min([Gap.CSDP2, Gap.CSDP3, Gap.CSDP4],[],2), Gap.PATH];
hb = bar(Y,'grouped');
ax2 = gca;
ax2.XTick = 1:size(TSPLIB,1);
ax2.YTick = 0:100:1000;
set(ax2,'xticklabel',TSPLIB(:,1))
ax2.TickLabelInterpreter = 'latex';
ax2.FontSize = 11;
xtickangle(45)
set(ax2, 'YGrid', 'on', 'XGrid', 'off')
set(gca,'LineWidth',0.5,'TickLength',[0 0]);
ylabel('relative gap \%','interpreter','latex','fontsize',13)

ax1.Position = ax2.Position;

%% Last edit: Alp Yurtsever - December 2, 2019
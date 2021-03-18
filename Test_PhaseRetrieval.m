%%  Test setup for Phase Retrieval SDP - by Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)

%% Choose data and method
% NOTE: Before the first time running this script, you should draw the
% random problem instances (create random data) by running "CreateData.m"
% in "FilesPhaseRetrieval" folder. 

L = 12;
METHOD = 'SketchyCGAL';
n = 1000;
MC = 1;
% METHOD = 'SketchyCGAL'; % 'CGAL', 'SketchyCGAL', 'ThinCGAL'
% n = 10000; % problem size, choose n from 1e2, 1e3, ..., 1e6
% MC = 1; % MonteCarlo iteration, choose MC from 1, 2, ... 20.

%% Preamble
rng(0,'twister');
addpath utils
addpath solver

%% Load data

load(['FilesPhaseRetrieval/',num2str(L),'/data/',num2str(n),'_',num2str(MC),'.mat']);

%% Construct the Phase Retrieval SDP

Aop_vec = @(x) evalAop_vec(x,Masks); % Forward operator (sampling operator)
AopTrPower = @(y,x) evalAopTr_pow(y,x,n,L,Masks); % Backword MULTIPLICATION operator
N = size(b,1); % total # of measurements

%% Shared parameters
% lambda0 = 1; % we didn't tune - choose 1 - you can tune this!
% maxit = 1e6; % limit on number of iterations
% saveit = round([2.^(0:log2(maxit)),maxit]); % error will be measured only at these iterations (and the last iteration)

nxt = norm(x_true); % norm of the solution
stopCond = @(x) norm(x_true - exp(-1i*angle(x_true'*x(:,1)))*x(:,1))/nxt; % stopping criterion (relative distance to solution)

%% Construct the black-box oracles

Primitive1 = @(x) x;
Primitive2 = @(y,x) evalAopTr_pow(y,x,n,L,Masks);
Primitive3 = @(x) evalAop_vec(x,Masks);
a = [0, 3*n];

% Compute scaling factors
SCALE_X = 1/max(a);
SCALE_C = 1/sqrt(n); % 1/norm(C,'fro')

Aoper = @(x) reshape(fft(bsxfun(@times,conj(Masks),x)),[],1)/sqrt(n);
AToper = @(x) sum(Masks.*(ifft(reshape(x,n,L))*sqrt(n)),2);
[x,~] = eigs(@(x) AToper(Aoper(x)), n, 1, 'LM');
NORM_A_ESTIMATE = norm(Aop_vec(x));
SCALE_A = 1/NORM_A_ESTIMATE; 
clearvars x

%% Debugging mode: Norm of A

% F = dftmtx(n)/sqrt(n);
% Adebug = nan(n^2,N);
% ptr = 0;
% for l = 1:size(Masks,2)
%     for t = 1:n
%         ptr = ptr+1;
%         f = Masks(:,l) .* conj(F(:,t)); 
%         Adebug(:,ptr) = reshape(f*f',[],1); 
%     end
% end
% Adebug = Adebug';
% 
% x = randn(n,1) + 1i*randn(n,1);
% x = x/norm(x);
% asd = (Adebug*vec(x*x') - Primitive3(x));
% max(abs(real(asd)))
% max(abs(imag(asd)))

%% Debugging mode: Check transpose

% xDebug = randn(n,1);
% yDebug = randn(size(b,1),1);
% ip1 = yDebug'*Primitive3(xDebug);
% ip2 = Primitive2(yDebug,xDebug)'*xDebug;
% if norm(ip1 - ip2) > 1e-10 % should be small
%     warning('Please check the black box oracles!');
% end

%% Start memory logging
% NOTE: This works only on Unix systems. 

hmL = memLog(['PhaseRetrieval_',num2str(L),'_',num2str(n),'_',num2str(MC),'_',METHOD]);
hmL.start;
MEMBEGIN = hmL.prompt;

%% Solve with SketchyCGAL
R = 5; % rank/sketch size parameter
beta0 = 1; % we didn't tune - choose 1 - you can tune this!
K = inf;
maxit = 1e6; % limit on number of iterations

timer = tic;
cputimeBegin = cputime;

out = CGAL( n, Primitive1, Primitive2, Primitive3, a, b, R, maxit, beta0, K, ...
    'field','complex',...
    'METHOD',METHOD,...
    'SCALE_C',SCALE_C,...
    'SCALE_X',SCALE_X,...
    'SCALE_A',SCALE_A,...
    'FLAG_MULTRANK_P1',true,... % err defines the spectral rounding for maxcut
    'FLAG_MULTRANK_P3',false,... % err defines the spectral rounding for maxcut
    'stopcond',stopCond,'stoptol',1e-2); % algorithm stops when 1e-2 relative distance to solution is reached

cputimeEnd = cputime;
totalTime = toc(timer);

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

%% Stop memory logging

MEMEND = hmL.prompt;
hmL.stop;
out.memory = (MEMEND - MEMBEGIN)/1000; %in MB

%% Save results

if ~exist(['./results/PhaseRetrieval/',num2str(L)],'dir'), mkdir(['./results/PhaseRetrieval/',num2str(L)]); end
save(['./results/PhaseRetrieval/',num2str(L),'/',num2str(n),'_',num2str(MC),'_',METHOD],'out','-v7.3');

%% Oracle Implementations
% NOTE: Defining a function in a MATLAB script was not available in older
% versions. If you are using an old version of MATLAB, you might need to
% save this function as a seperate ".m" file in your path.

function out = evalAop_vec(x,Masks)

Ax = bsxfun(@times,conj(Masks),x);
Ax = fft(Ax)/sqrt(size(Ax,1));
Ax = abs(Ax).^2;
out = Ax(:);

end

function out = evalAopTr_pow(y,x,n,L,Masks)

Ax = bsxfun(@times,conj(Masks),x);
Ax = fft(Ax)/sqrt(n);

Atop = reshape(y,n,L) .* Ax;
Atop = ifft(Atop)*sqrt(n);
Atop = Masks.*Atop;
out = sum(Atop,2);

end

%% Last edit: Alp Yurtsever - July 24, 2019
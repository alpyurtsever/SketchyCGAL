%% Test setup for Ptychography SDP - by Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)
% NOTE: You need to run "./FilesPtychography/CreateData.m" to generate
% draw a rand transmission matrix for the ptychogrpahy system.

%% Choose data
dim = 320; % choose 160, 320 or 640

%% Preamble
rng(0,'twister');
addpath FilesPtychography
addpath solver;

%% Load data
% NOTE: We used the transmission matrix model of a working ptychography 
% system in our experiments. However, due to the restrictions on the data, 
% this script uses a synthetic model instead. Hence, the results will be 
% different than the ones in Figure 7.4 and Figure SM6.2 in our paper. 
% We generate the synthetic model by randomly drawing the complex 
% coefficients of the transmission matrix (with unit magnitude and
% i.i.d. random phase angle chosen uniformly over -pi to pi.) 
load(['FilesPtychography/data/SyntheticTrMat',num2str(dim),'.mat']);

n = size(A,2);
N = size(A,1);

vec         = @(x) x(:);
F           = @(x) vec(ifft2(ifftshift(reshape(x, mask_width,mask_width,num_im)))*mask_width);
FT          = @(x) vec(fftshift(fft2(reshape(x, mask_width,mask_width,num_im))/mask_width));
Aoper       = @(x) F(A*x);
AToper      = @(x) A'*FT(x);
Aop_vec     = @(x) abs(Aoper(x)).^2;
AopTrPower  = @(y,x) AToper(y(:) .* Aoper(x));

%% Read image data
imtrue = imread(['FilesPtychography/data/Cell',num2str(dim),'.tiff']); % original image
imtrueFourier = fftshift(fft2(imtrue)); % original image in Fourier domain
imtrueFourier = imtrueFourier./norm(imtrueFourier,'fro'); % normalization

%% Construct the Ptychography SDP

Primitive1 = @(x) x;
Primitive2 = @(y,x) AToper(y(:) .* Aoper(x));
Primitive3 = @(x) abs(Aoper(x)).^2;
a = [0, 1.5];

% Compute scaling factors
SCALE_X = 1/max(a);
SCALE_C = 1/sqrt(n); % 1/norm(C,'fro')

[x,~] = eigs(@(x) AToper(Aoper(x)), n, 1, 'LM');
NORM_A_ESTIMATE = norm(Aop_vec(x));
SCALE_A = 1/NORM_A_ESTIMATE; 
clearvars x
% NOTE: This is a very rough estimation of the order or the norm of A. Note
% that we computed the eigenvector in the n-dimensional vector space, not
% in the (n x n)-dimensional space. I am not sure how we can get a good
% approximation of norm of A without using O(n^2) storage. If we could
% formulate A_i as sparse matrices, then we could compute the norm of A*A'
% by using O(d) storage, if A*A' is also sparse enough. In this case, A_i
% are easy to compute with fft, but we cannot generate sparse A_i. In any
% case, I hope that a rough estimate to work well enough for pre-scaling.
% One can tune beta0 for better performance.

%% Simulate FP measurements
b = Primitive3(imtrueFourier(:));

%% Start memory logging
% NOTE: This works only on Unix systems. 

hmL = memLog(['Ptychography_',num2str(dim)]);
hmL.start;
MEMBEGIN = hmL.prompt;

%% Solve

beta0 = 1; % we didn't tune - choose 1 - you can tune this!
K = inf;
R = 5;
maxit = 1e4; % limit on number of iterations

% saveit = [10,100,1000,10000]; % error will be measured only at these iterations (and the last iteration)
saveit = [1,2,4,6,8,10:10:90,100:100:900,1000:1000:10000];

err = {};
err{1} = 'Plotter';
err{2} = @(x) RoundAndPlot(x, dim, saveit);

timer = tic;
cputimeBegin = cputime;

out = CGAL( n, Primitive1, Primitive2, Primitive3, a, b, R, maxit, beta0, K, ...
    'field','complex',...
    'SCALE_C',SCALE_C,...
    'SCALE_A',SCALE_A,...
    'SCALE_X',SCALE_X,...
    'FLAG_MULTRANK_P1',true,... % err defines the spectral rounding for maxcut
    'FLAG_MULTRANK_P3',false,... % err defines the spectral rounding for maxcut
    'errfncs',err,... % err defines the spectral rounding for maxcut
    'saveit',saveit); % we collect data about the convergence at iterations 1,2,4,8,...,T

cputimeEnd = cputime;
totalTime = toc(timer);

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

%% Stop memory logging

MEMEND = hmL.prompt;
hmL.stop;
out.memory = (MEMEND - MEMBEGIN)/1000; %in MB

%% Save results

if ~exist('./results/Ptychography','dir'), mkdir('./results/Ptychography'); end
save(['./results/Ptychography/',num2str(dim),'_info.mat'],'out','-v7.3');

%% Implement rounding
% NOTE: Defining a function in a MATLAB script was not available in older
% versions. If you are using an old version of MATLAB, you might need to
% save this function as a seperate ".m" file in your path.

function out = RoundAndPlot( xestimate, dim, saveit )
% Implements a spectral rounding and shows the constructed image. 

persistent ptr
if isempty(ptr)
    ptr = 0;
end
ptr = ptr+1;

n = sqrt(size(xestimate,1));

Img = abs(ifft2(ifftshift(ifftshift(reshape(xestimate(:,1),n,n),1),2)));
Img = Img - min(Img(:));
Img = Img./max(Img(:));
Img = uint8(255*abs(Img));

if ~exist('./results/Ptychography','dir'), mkdir('./results/Ptychography'); end
% imwrite(Img, ['./results/Ptychography/',num2str(dim),'_',num2str(saveit(ptr)),'.tiff']);
% imwrite(Img, ['./results/Ptychography/',num2str(dim),'_',num2str(saveit(ptr)),'.png']);
imshow(Img); drawnow;

out = nan;

end

%% Last edit: Alp Yurtsever - December 06, 2019

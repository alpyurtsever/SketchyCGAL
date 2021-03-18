%% Test setup for QAP SDP - by Alp Yurtsever (alp.yurtsever@epfl.ch - alpy@mit.edu)

%% Choose data
% NOTE: You need to download data from QAPLIB and locate them to under the
% "FilesQAP/data/qapdata/" folder (resp. TSPLIB, "FilesQAP/data/tspdata/").
% Links for QAPLIB: 
%   http://anjos.mgi.polymtl.ca/qaplib/inst.html
%   https://www.opt.math.tugraz.at/qaplib/inst.html
% Links for TSPLIB: 
%   http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/index.html
%   http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95

% dataname = 'qapdata/chr12a'; % you can choose other data files as well
% dataname = 'tspdata/ulysses16';

%% Preamble
rng(0,'twister');
addpath solver
addpath utils;
addpath FilesQAP;
addpath FilesQAP/munkres;
addpath FilesQAP/PartialTrace;

%% Load data

if strncmp(dataname,'qap',3)
    [A,B] = qapread(['./FilesQAP/',dataname,'.dat']);
elseif strncmp(dataname,'tsp',3)
    [A,B] = tspread(['./FilesQAP/',dataname,'.tsp']);
else
    error('Unknown data type.')
end

k = size(A,1);
n = k^2 + 1;

if nnz(A)/numel(A) <= 0.25
    A = sparse(A);
end

if nnz(B)/numel(B) <= 0.25
    B = sparse(B);
end

if nnz(A) <= nnz(B)
    AA = A;
    A = B;
    B = AA;
    clearvars AA;
end % As convention we pick the sparse one as B 


%% Useful functionals for efficient implementation of Primitives
% Note that X = [1, vec(P)'; vec(P), Y];

% vectorization
vec = @(x) x(:);
% y = AkronBx(A,B,x) -> Computes kron(A,B)*x
AkronBx = @(A,B,x) vec(B*reshape(x,size(B,2),size(A,2))*A');
% y = Pvec(x) -> Samples P from decision variable and vectorizes it
Pvec = @(x) x(1).*x(2:end);
% y = P(x) -> Samples P from decision variable
P = @(x) reshape(Pvec(x),[k,k]);
% y = diagY(x) -> Samples diagonal of Y from decision variable
diagY = @(x) x(2:end).^2;
% y = Yvec(x) -> Samples Y from decision variable and vectorizes it
Yvec = @(x) vec(x(2:end)*x(2:end)');

% find indexes where nonzeros should be enforced (we could implement this more efficiently!)
Bpattern = (B ~= 0 + speye(size(B))) ~= 0;
Cpattern = kron(Bpattern,ones(k,k));
Cpattern = [0, zeros(1,size(Cpattern,2)); zeros(size(Cpattern,1),1), Cpattern];
indNNZ = find(Cpattern);
clearvars Bpattern Cpattern;
nNNZ = numel(indNNZ);
[iNNZ,jNNZ] = ind2sub([n,n],indNNZ);

%% Construct QAP SDP

% Constraint: X(1,1,) = 1
% Constraint: Diag(Y) = vec(P) 
% Constraint: P1 = 1 
% Constraint: 1'P = 1'
% Trace1 constraint
% Trace2 constraint
% Nonnegativity constraints

% We use "TrX" function by Tony Cubitt to implement Trace1 and Trace2 constraints. 
% See "FilesQAP/PartialTrace/TrX.m" for the code and the disclaimer.
% Link: http://www.dr-qubit.org/matlab.html

%% Approach 1

Primitive1 = @(x) [0;AkronBx(B,A,x(2:end))];

ind1 = 1;
ind2 = 2:(k^2+1);
ind3 = (k^2+1+1):(k^2+1+k);
ind4 = (k^2+1+k+1):(k^2+1+k+k);
ind5 = (k^2+1+k+k+1):(k^2+1+k+k+k^2);
ind6 = (k^2+1+k+k+k^2+1):(k^2+1+k+k+k^2+k^2);
ind7 = (k^2+1+k+k+k^2+k^2+1):(k^2+1+k+k+k^2+k^2+nNNZ);

% Matlab implementation
Primitive2_MATLAB = @(y,x) ...
    [y(ind1)*x(1); zeros(n-1,1)] + ...
    ([0;y(ind2)].*x - 0.5*([0;y(ind2)]*x(1) + [[0,y(ind2)']*x;zeros(n-1,1)])) + ...
    0.5*([0;repmat(y(ind3)*x(1),[k,1])] + [repmat(y(ind3),[k,1])'*x(2:end);zeros(k^2,1)]) + ...
    0.5*([0;repelem(y(ind4)*x(1),k)] + [repelem(y(ind4),k)'*x(2:end);zeros(k^2,1)]) + ...
    [0; vec(reshape(y(ind5),[k,k])*reshape(x(2:end),[k,k]))] + ...
    [0; vec(reshape(x(2:end),k,k)*reshape(y(ind6),[k,k])')] + ...
    mexSparseMult(iNNZ,jNNZ,y(ind7),n,x);

% Matlab implementation
Primitive3_MATLAB = @(x) [...
    x(1)^2; ...
    (diagY(x) - Pvec(x)); ...
    sum(P(x),2); ...
    sum(P(x),1)'; ...
    vec(TrX(x(2:end),1,[k,k]));...
    vec(TrX(x(2:end),2,[k,k]));...
    x(iNNZ).*x(jNNZ) ];

a = k+1;

b1 = 1;
b2 = zeros(k^2,1);
b3 = ones(k,1);
b4 = ones(k,1);
b5 = vec(eye(k));
b6 = vec(eye(k));
b7 = [zeros(nNNZ,1),ones(nNNZ,1)];

b = [b1,b1;...
    b2,b2;...
    b3,b3;...
    b4,b4;...
    b5,b5;...
    b6,b6;...
    b7];

clearvars b1 b2 b3 b4 b5 b6 b7

%% Approach 2 - with MEX

% Constraint: X(1,1,) = 1
subN = 1;
subI = 1;
subJ = 1;
VAL = 1;

% Diag(Y) = vec(P) constraint
J = ones(2*(n-1),1);
J(2:2:2*n-1) = (2:n)';
I = repelem((2:n)', 2, 1);
N = repelem((1:n-1)', 2, 1);
V = ones(2*(n-1),1);
V(1:2:end) = -1;

subJ = [subJ; J];
subI = [subI; I];
subN = [subN; N + subN(end)];
VAL = [VAL; V];

% P1 = 1 constraint
I = ones(n-1,1);
J = reshape(reshape(2:n,k,k)',[],1);
N = repelem((1:k)', k, 1);
V = ones(n-1,1);

subJ = [subJ; J];
subI = [subI; I];
subN = [subN; N + subN(end)];
VAL = [VAL; V];

% 1'P = 1' constraint
I = ones(n-1,1);
J = (2:n)';
N = repelem((1:k)', k, 1);
V = ones(n-1,1);

subJ = [subJ; J];
subI = [subI; I];
subN = [subN; N + subN(end)];
VAL = [VAL; V];

% Trace1 constraint
I = reshape(repelem(reshape((2:n),k,k)',1,k),[],1);
J = vec(repmat(reshape((2:n),k,k)',1,k));
N = repelem((1:n-1)', k, 1);
V = ones(k^3,1);

subJ = [subJ; J];
subI = [subI; I];
subN = [subN; N + subN(end)];
VAL = [VAL; V];

% % Trace2 constraint
I = reshape(repelem(reshape((2:n),k,k),1,k),[],1);
J = vec(repmat(reshape((2:n),k,k),1,k));
N = repelem((1:n-1)', k, 1);
V = ones(k^3,1);

subJ = [subJ; J];
subI = [subI; I];
subN = [subN; N + subN(end)];
VAL = [VAL; V];

% % Nennegativity constraints
subI = [subI; iNNZ];
subJ = [subJ; jNNZ];
subN = [subN; (1:nNNZ)' + subN(end)];
VAL = [VAL; ones(nNNZ,1)];

uindN = int32(subN);
uindI = int32(subI);
uindJ = int32(subJ);
m = double(max(subN));

Primitive2_MEX = @(y,x) mexPrimitive2(uindN,uindI,uindJ,VAL,y,x);
Primitive3_MEX = @(x) mexPrimitive3(uindN,uindI,uindJ,VAL,m,x);

%% Test to find the faster approach
x = randn(n,1);

% Check for Primitive 3
tic;
for t = 1:30
    y = Primitive3_MEX(x);
end
MEXTIME = toc;
tic;
for t = 1:30
    yy = Primitive3_MATLAB(x);
end
if norm(y-yy) > 1e-9
    error('Primitive3_MEX and Primitive3_MATLAB results are inconsistent.')
end
MATLABTIME = toc;
if MEXTIME < MATLABTIME
    Primitive3 = Primitive3_MEX;
    fprintf('MEX implementation of Primitive 3 is chosen.\n')
else
    Primitive3 = Primitive3_MATLAB;
    fprintf('MATLAB implementation of Primitive 3 is chosen.\n')
end
clearvars Primitive3_MEX Primitive3_MATLAB

% Check for Primitive 2
tic;
for t = 1:30
    xb = Primitive2_MEX(y,x);
end
MEXTIME = toc;
tic;
for t = 1:30
    xxb = Primitive2_MATLAB(y,x);
end
if norm(xb-xxb) > 1e-9
    error('Primitive2_MEX and Primitive2_MATLAB results are inconsistent.')
end
MATLABTIME = toc;
if MEXTIME < MATLABTIME
    Primitive2 = Primitive2_MEX;
    fprintf('MEX implementation of Primitive 2 is chosen.\n')
else
    Primitive2 = Primitive2_MATLAB;
    fprintf('MATLAB implementation of Primitive 2 is chosen.\n')
end
clearvars Primitive3_MEX Primitive3_MATLAB

%% Find the SCALE factors

SCALE_X = 1/(k+1);
SCALE_C = 1/sqrt((norm(A,'fro')*norm(B,'fro'))^2+1);

IND = sub2ind([n,n],subI,subJ);
Amatrix = sparse(subN, IND, VAL, m, n^2); % Vectorized matrix representation of A --> b = Amatrix * vec(X)
SCALE_A = 1 ./ norms(Amatrix,2,2);

Amatrix = spdiags(SCALE_A,0,m,m)*Amatrix;
Amatrix = Amatrix*Amatrix';
normA = sqrt(eigs(Amatrix,1,'LM'));
clearvars Amatrix;
SCALE_A = SCALE_A ./ normA;

%% Implement rounding

err = {};
err{1} = 'eigsRound'; % name for error
err{2} = @(u) eigsRound(A,B,u); % function definition at the bottom of this script

%% Start memory logging
% NOTE: This works only on Unix systems. 

hmL = memLog(dataname);
hmL.start;
MEMBEGIN = hmL.prompt;

%% Run SketchyCGAL

beta0 = 1; % we didn't tune - choose 1 - you can tune this!
K = inf;
maxit = 1e6; % limit on number of iterations
R = k; % rank/sketch size parameter

timer = tic;
cputimeBegin = cputime;

out = CGAL( n, Primitive1, Primitive2, Primitive3, a, b, R, maxit, beta0, K, ...
    'SCALE_A',SCALE_A,...
    'SCALE_C',SCALE_C,...
    'SCALE_X',SCALE_X,...
    'errfncs',err,... % err defines the spectral rounding for maxcut
    'walltime',3*24*60*60); % algorithm stops after at most 3 days

totalTime = toc(timer);
cputimeEnd = cputime;

out.totalTime = totalTime;
out.totalCpuTime = cputimeEnd - cputimeBegin;

%% Stop memory logging

MEMEND = hmL.prompt;
hmL.stop;
out.memory = (MEMEND - MEMBEGIN)/1000; %in MB

%% Save Results

if ~exist(['results/QAP/',dataname],'dir'), mkdir(['results/QAP/',dataname]); end
save(['results/QAP/',dataname,'/SketchyCGAL.mat'],'out','-v7.3');

%% Implement rounding
% NOTE: Defining a function in a MATLAB script was not available in older
% versions. If you are using an old version of MATLAB, you might need to
% save this function as a seperate ".m" file in your path.

% We use the MUNKRES (aka Hungarian Algorithm) implementation by Yi Ciao. 
% See "FilesQAP/munkres" folder for the code and the disclaimer. 
% Link: https://www.mathworks.com/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems-v2-3

function out = eigsRound(A,B,x)
out = inf;
k = size(A,1);
for t=1:size(x,2)
    P = reshape(x(2:end,t),[k,k]);
    order = munkres(max(P(:))-P);
    P = sparse(1:k,order,ones(k,1));
    out = min(out,trace(A*P*B*P'));
    % Also try minus P (note that the sign of the eigenvector is ambigous!)
    P = reshape(-x(2:end,t),[k,k]);
    order = munkres(max(P(:))-P);
    P = sparse(1:k,order,ones(k,1));
    out = min(out,trace(A*P*B*P'));
end
end

%% Last edit: Alp Yurtsever - July 24, 2020
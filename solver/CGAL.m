function [ out, U, Delt ] = CGAL( n, Primitive1, Primitive2, Primitive3, a, b, R, T, beta0, K, varargin )
%CGAL This function implements CGAL - SketchyCGAL - ThinCGAL for solving
%semidefinite programs of the following form:
%
%      minimize  <C,X>  subj.to   X is symmetric positive semidefinite
%                                 al <= tr(X) <= au
%                                 bl_i <= <A_i,X> <= bu_i
%
% Coded by: Alp Yurtsever - alp.yurtsever@epfl.ch
% Last modified: 19 November, 2019

narginchk(8,inf);
if isempty(beta0), beta0 = 1; end
if isempty(K), K = inf; end

% User defined options!
WALLTIME = inf; % Default: Walltime is infinite
STOPTOL = []; STOPCOND = []; % Default: No stopping condition
FLAG_INCLUSION = ~isvector(b); % [true, false] = ['AX = b', 'b(:,1) <= AX <= b(:,2)]
FLAG_LANCZOS = true; % [true, false] = [Lanczos Method, Power Method]
FLAG_TRACECORRECTION = true; % optional trace correction step
FLAG_METHOD = 1; % [0,1,2] = [CGAL, SketchyCGAL, ThinCGAL]
FLAG_MULTRANK_P1 = false;
FLAG_MULTRANK_P3 = false;
SKETCH_FIELD = 'real'; % 'real' or 'complex'
ERR = {};
SAVEHIST = unique([2.^(0:floor(log2(T))),T])';

SCALE_A = 1; % can be a scale or vector of size b
SCALE_C = 1; % a scale only
SCALE_X = 1; % a scale only
NORM_A = 1; % Ideally, this should be norm of A

if ~isempty(varargin)
    for tt = 1:2:length(varargin)
        switch lower(varargin{tt})
            case 'method'
                switch lower(varargin{tt+1})
                    case 'cgal'
                        FLAG_METHOD = 0;
                    case 'sketchycgal'
                        FLAG_METHOD = 1;
                    case 'thincgal'
                        FLAG_METHOD = 2;
                    otherwise
                        error(['Unknown optimization method: ', varargin{tt+1}]);
                end
            case 'walltime'
                WALLTIME = varargin{tt+1};
            case 'errfncs'
                ERR = varargin{tt+1};
            case 'flag_multrank_p1'
                FLAG_MULTRANK_P1 = logical(varargin{tt+1});
            case 'flag_multrank_p3'
                FLAG_MULTRANK_P3 = logical(varargin{tt+1});
            case 'norma'
                NORM_A = varargin{tt+1};
            case 'scale_a'
                SCALE_A = varargin{tt+1};
            case 'scale_c'
                SCALE_C = varargin{tt+1};
            case 'scale_x'
                SCALE_X = varargin{tt+1};
            case 'stoptol'
                STOPTOL = varargin{tt+1}; % NOTE: Only for "AX = b" type constraints for now!!!
            case 'stopcond'
                STOPCOND = varargin{tt+1};
            case 'saveit'
                SAVEHIST = varargin{tt+1};
            case 'field'
                SKETCH_FIELD = varargin{tt+1};
            case 'tracecorrection'
                FLAG_TRACECORRECTION = logical(varargin{tt+1});
            case 'lmo'
                switch lower(varargin{tt+1})
                    case 'power'
                        FLAG_LANCZOS = false;
                    case 'lanczos'
                        FLAG_LANCZOS = true;
                    otherwise
                        error('Unknown linear minimization oracle type.');
                end
            otherwise
                warning(['Unknown option: ',varargin{tt}]);
        end
    end
end

% Scale the problem
b_org = b;
a_org = a;
RESCALE_OBJ = 1;
RESCALE_FEAS = 1;
if any(SCALE_A ~= 1)
    b = b .* SCALE_A;
    Primitive3 = @(x) Primitive3(x) .* SCALE_A;
    Primitive2 = @(y,x) Primitive2(y .* SCALE_A, x);
    RESCALE_FEAS = RESCALE_FEAS ./ SCALE_A;
end
if SCALE_X ~= 1
    b = b .* SCALE_X;
    a = a .* SCALE_X;
    RESCALE_OBJ = RESCALE_OBJ / SCALE_X;
    RESCALE_FEAS = RESCALE_FEAS / SCALE_X;
end
if SCALE_C ~= 1
    Primitive1 = @(x) Primitive1(x) * SCALE_C;
    RESCALE_OBJ = RESCALE_OBJ / SCALE_C;
end

if FLAG_INCLUSION
    PROJBOX = @(y) min(max(y,b(:,1)),b(:,2));
end
% projK = @(y) projK(y.*scaleFeas)./scaleFeas; %check this

% Create OUT stuct where we store optimization information
[out, errNames, errNamesPrint, ptr] = createErrStructs();


% Initialize the decision variable and the dual
if FLAG_METHOD == 0
    X = zeros(n,n);
elseif FLAG_METHOD == 1
    mySketch = NystromSketch(n, R, SKETCH_FIELD);
elseif FLAG_METHOD == 2
    UTHIN = zeros(n,1);
    DTHIN = 0;
else
    error('Unknown approach for storing decision variable.');
end

% Initialize the dual
z = zeros(size(b,1),1);
y0 = zeros(size(b,1),1);
y = y0;
objective = 0;

% Choose the linear minimization oracle
if FLAG_LANCZOS
    ApproxMinEvec = @(x,t) ApproxMinEvecLanczos(x, n, ceil((t^0.25)*log(n))); % More practical (Almost-storage optimal but faster)
else
    ApproxMinEvec = @(x,t) ApproxMinEvecPower(x, n, ceil(8*(t^0.5)*log(n))); % Storage-optimal but slower
end

cntTotal = 0; % Oracle counter for Primitive2

% Start the timer
timer = tic;
cputime0 = cputime;
totTime = toc(timer);
totCpuTime = cputime - cputime0;

TRACE = 0;

for t = 1:T
    
    beta = beta0*sqrt(t+1);
    eta = 2/(t+1);
    
    if FLAG_INCLUSION
        vt = y + beta.*(z - PROJBOX(z + (1/beta).*y));
    else
        vt = y + beta*(z-b);
    end
    eigsArg = @(u) Primitive1(u) + Primitive2(vt,u);
    
    [u,sig,cntInner] = ApproxMinEvec(eigsArg,t);
    cntTotal = cntTotal + cntInner;
    
    if sig > 0, a_t = min(a); else, a_t = max(a); end
    u = sqrt(a_t)*u;
    
    TRACE = (1-eta)*TRACE + eta*a_t;
    
    % Check stopping criteria
    if ~isempty(STOPTOL)
        if ~isempty(STOPCOND)
            switch FLAG_METHOD
                case 0
                    [U,Delt] = eigs(X, R, 'LM');
                case 1
                    [U,Delt] = mySketch.Reconstruct();
                    if FLAG_TRACECORRECTION
                        Delt = Delt + ((TRACE-trace(Delt))/R)*eye(size(Delt));
                    end
                case 2
                    U = UTHIN;
                    Delt = DTHIN;
            end
            U = U*sqrt(Delt);
            if STOPCOND(U / sqrt(SCALE_X)) <= STOPTOL
                if ~isempty(SAVEHIST), updateErrStructs(); printError(); clearErrStructs(); end
                out.status = 'stopping criteria met (stopCond)';
                break;
            end
        else
            FeasCond = norm(z - b);
            ObjCond = objective + (y + beta.*(z-b))'*b + 0.5*beta*FeasCond^2 - eigsArg(u)'*u;
            if FeasCond <= STOPTOL && ObjCond <= STOPTOL
                if ~isempty(SAVEHIST), updateErrStructs(); printError(); clearErrStructs(); end
                out.status = 'stopping criteria met';
                break;
            end
        end
    end
    
    zEvec = Primitive3(u);
    z = (1-eta)*z + eta*zEvec;
    
    objEvec = u'*Primitive1(u);
    objective = (1-eta)*objective + eta*objEvec;
    
    if FLAG_METHOD == 0
        X = (1-eta)*X + eta*(u*u');
    elseif FLAG_METHOD == 1
        mySketch.RankOneUpdate(u,eta);
    elseif FLAG_METHOD == 2
        [UTHIN, DTHIN] = updateThinSVDrank1(UTHIN,(1-eta)*DTHIN,u,eta);
    else
        error('Unknown approach for storing decision variable.');
    end
    
    beta_p = beta0*sqrt(t+2);
    
    if FLAG_INCLUSION
        dualUpdate = z - PROJBOX( z + (1/beta_p).*y );
    else
        dualUpdate = z - b;
    end
    
    sigma = min(beta0, 4 * beta_p * eta^2 * max(a)^2 * NORM_A^2 / norm(dualUpdate)^2);
    
    % Update the DUAL
    yt1 = y + sigma.*dualUpdate;
    if norm(yt1 - y0) <= K, y = yt1; end
    
    % Measure the runtime
    totTime = toc(timer);
    totCpuTime = cputime - cputime0;
    if totCpuTime > WALLTIME
        if ~isempty(SAVEHIST), updateErrStructs(); printError(); clearErrStructs(); end
        out.status = 'wall time achieved';
        break;
    end
    
    % Update OUT
    if any(t==SAVEHIST), updateErrStructs(); printError(); end
    
end

if t == T
    out.status = 'maximum number of iterations achieved';
end

if FLAG_METHOD == 0
    U = X / SCALE_X;
    Delt = [];
elseif FLAG_METHOD == 1
    if FLAG_TRACECORRECTION
        Delt = Delt + ((TRACE-trace(Delt))/R)*eye(size(Delt));
    end
    Delt = Delt / SCALE_X;
elseif FLAG_METHOD == 2
    U = UTHIN;
    Delt = DTHIN / SCALE_X;
end


%% Nested functions
    function [out, errNames, errNamesPrint, ptr] = createErrStructs()
        if ~isempty(STOPTOL)
            if ~isempty(STOPCOND)
                out.info.stopCond = nan(numel(SAVEHIST),1);
            else
                out.info.stopObj = nan(numel(SAVEHIST),1);
                out.info.stopFeas = nan(numel(SAVEHIST),1);
            end
        end
        out.info.primalObj = nan(numel(SAVEHIST),1);
        out.info.primalFeas = nan(numel(SAVEHIST),1);
        if FLAG_METHOD == 1
            out.info.skPrimalObj = nan(numel(SAVEHIST),1);
            out.info.skPrimalFeas = nan(numel(SAVEHIST),1);
        end
        for eIt = 1:2:length(ERR)
            out.info.(ERR{eIt}) = nan(numel(SAVEHIST),1);
        end
        out.info.cntTotal = nan(numel(SAVEHIST),1);
        out.iteration = nan(numel(SAVEHIST),1);
        out.time = nan(numel(SAVEHIST),1);
        out.cputime = nan(numel(SAVEHIST),1);
        errNames = fieldnames(out.info);
        
        errNamesPrint = errNames;
        for pIt = 1:length(errNamesPrint)
            if length(errNamesPrint{pIt}) > 10
                errNamesPrint{pIt} = errNamesPrint{pIt}(1:10);
            end
        end
        ptr = 0;
        
        out.params.ALPHA = a;
        out.params.BETA0 = beta0;
        out.params.R = R;
        out.params.FLAG_LANCZOS = FLAG_LANCZOS;
        out.params.FLAG_TRACECORRECTION = FLAG_TRACECORRECTION;
        out.params.FLAG_SKETCH = FLAG_METHOD;
        out.status = 'running';
    end

    function updateErrStructs()
        if mod(ptr,20) == 0, printHeader(); end
        ptr = ptr+1;
        out.iteration(ptr,1) = t;
        out.time(ptr,1) = totTime;
        out.cputime(ptr,1) = totCpuTime;
        out.info.primalObj(ptr,1) = objective * RESCALE_OBJ;
        if FLAG_INCLUSION, FEAS = norm((z - PROJBOX(z)) .* RESCALE_FEAS); % maybe / (norm(PROJBOX(z) .* RESCALE_FEAS) + 1); ?
        else, FEAS = norm((z - b) .* RESCALE_FEAS) / max(norm(b_org),1);
        end
        out.info.primalFeas(ptr,1) = FEAS;
        out.info.cntTotal(ptr,1) = cntTotal;
        if FLAG_METHOD == 1
            [U,Delt] = mySketch.Reconstruct();
            if FLAG_TRACECORRECTION
                Delt = Delt + ((TRACE-trace(Delt))/R)*eye(size(Delt));
            end
            U = U*sqrt(Delt);
            out.info.skPrimalObj(ptr,1) = trace(U'*Primitive1MultRank(U)) * RESCALE_OBJ;
            if FLAG_INCLUSION
                AUU = Primitive3MultRank(U);
                out.info.skPrimalFeas(ptr,1) = norm((AUU - PROJBOX(AUU)) .* RESCALE_FEAS); % maybe /  (norm(PROJBOX(AUU) .* RESCALE_FEAS) + 1); ?
            else
                out.info.skPrimalFeas(ptr,1) = norm((Primitive3MultRank(U) - b) .* RESCALE_FEAS) / max(norm(b_org),1);
            end
                %norm((Ax_b_sketch - projK(Ax_b_sketch)).*scaleFeas);
        elseif FLAG_METHOD == 0
            [U,Delt] = eigs(X, R, 'LM');
            U = U*sqrt(Delt);
        elseif FLAG_METHOD == 2
            [~,inds] = sort(diag(DTHIN),'descend');
            RR = min(R,size(UTHIN,2));
            inds = inds(1:RR);
            U = UTHIN(:,inds);
            Delt = DTHIN(1:RR,1:RR);
            U = U*sqrt(Delt);
        else
            error('Unknown FLAG_METHOD.');
        end
        if ~isempty(STOPTOL)
            if ~isempty(STOPCOND)
                out.info.stopCond(ptr,1) = STOPCOND(U / sqrt(SCALE_X)); %% alphaorg
            else
                out.info.stopObj(ptr,1) = ObjCond;
                out.info.stopFeas(ptr,1) = FeasCond;
            end
        end
        for eIt = 1:2:length(ERR)
            if nargin(ERR{eIt+1}) == 1
                out.info.(ERR{eIt})(ptr,:) = ERR{eIt+1}(U / sqrt(SCALE_X)); %% alphaorg
            end
        end
    end

    function clearErrStructs()
        out.iteration(isnan(out.iteration)) = [];
        out.time(isnan(out.time)) = [];
        out.cputime(isnan(out.cputime)) = [];
        fields = fieldnames(out.info);
        for i = 1:numel(fields)
            out.info.(fields{i})(isnan(out.info.(fields{i}))) = [];
        end
    end

    function printError()
        fprintf('| %7d |', t );
        for pIt = 1:length(errNames)
            if ~isnan(out.info.(errNames{pIt})(ptr))
                fprintf(' % 5.3e |', out.info.(errNames{pIt})(ptr) );
            else
                fprintf(' % 10s |', 'NaN' );
            end
        end
        fprintf('\n');
    end

    function printHeader()
        if ptr == 0
            if FLAG_METHOD == 0, fprintf('- CGAL SDP Solver - Beta.V.0.0 \n');
            elseif FLAG_METHOD == 1, fprintf('- SketchyCGAL SDP Solver - Beta.V.0.0 \n');
            elseif FLAG_METHOD == 2, fprintf('- ThinCGAL SDP Solver - Beta.V.0.0 \n');
            else, error('Unknown approach for storing decision variable.');
            end
        end
        fprintf(repmat('-',11 + length(errNamesPrint)*13,1));
        fprintf('\n');
        fprintf('|    iter | ');
        for pIt = 1:length(errNamesPrint)
            fprintf('%10s | ', errNamesPrint{pIt});
        end
        fprintf('\n');
        fprintf(repmat('-',11 + length(errNamesPrint)*13,1));
        fprintf('\n');
    end

    function CU = Primitive1MultRank(U)
        if FLAG_MULTRANK_P1
            CU = Primitive1(U);
        else
            CU = nan(size(U));
            for ind = 1:size(U,2)
                CU(:,ind) = Primitive1(U(:,ind));
            end
        end
    end

    function AUU = Primitive3MultRank(U)
        if FLAG_MULTRANK_P3
            AUU = Primitive3(U);
        else
            AUU = zeros(size(b));
            for ind = 1:size(U,2)
                AUU = AUU + Primitive3(U(:,ind));
            end
        end
    end

end



%% Shifted power method
function [v,xi,cnt] = ApproxMinEvecPower( M, n, q )

% [xi, v] = ApproxMinEig(M, q)
% Uses q iterations of shifted power iteration to find a unit vector v that
% approximately minimizes the quadratic form defined by the Hermitian matrix M.
% The output satisfies xi = v' * M * v with norm(v) = 1
% q = 1 + log(n) / eps iterations suffices to achieve eps error (on average over random start)

% Add consistency checks here
if isnumeric(M)
    M = @(x) M*x;
end

% Coarse spectral norm estimate via power method
[Mnorm, cnt] = powermethod(M,n,0.1);
shift = 2*Mnorm;

% Random initialization
v = randn( n, 1 );
v = v / norm( v );

% Main algorithm
Mv = M(v);
cnt = cnt + 1;

for iter = 1 : q  % q iterations of power method with matrix (shift * I - M)
    
    v = shift * v - Mv;
    v = v / norm( v );
    
    Mv = M(v);
    cnt = cnt + 1;
    
end

xi = v' * Mv;

end

function [e, cnt] = powermethod(S,n,tol)
cnt = 0;
x = randn(n,1);
x = x/norm(x);
e = 1;
e0 = 0;
while abs(e-e0) > tol*e
    e0 = e;
    Sx = S(x);
    if nnz(Sx) == 0, Sx = rand(size(Sx),class(Sx)); end
    x = S(Sx);
    normx = norm(x);
    e = normx/norm(Sx);
    x = x/normx;
    cnt = cnt+1;
    if cnt > 100, warning('Power method did not converge'); break; end
end
end

%% Lanczos method
function [v, xi, i] = ApproxMinEvecLanczos(M, n, q)
% Approximate minimum eigenvector
% Vanilla Lanczos method

q = min(q, n-1);                    % Iterations < dimension!

if isnumeric(M), M = @(x) M*x; end

Q = zeros(n, q+1);                  % Lanczos vectors

aleph = zeros(q,1);                 % Diagonal Lanczos coefs
beth = zeros(q,1);                  % Off-diagonal Lanczos coefs

Q(:,1) = randn(n, 1);               % First Lanczos vector is random
Q(:,1) = Q(:,1) / norm(Q(:,1));

for i = 1 : q
    Q(:, i+1) = M ( Q(:, i) );				% Apply M to previous Lanczos vector
    aleph(i) = real(Q(:, i)' * Q(:, i+1));		% Compute diagonal coefficients
    
    if (i == 1)                     % Lanczos iteration
        Q(:, i+1) = Q(:, i+1) - aleph(i) * Q(:, i);
    else
        Q(:, i+1) = Q(:, i+1) - aleph(i) * Q(:, i) - beth(i-1) * Q(:, i-1);
    end
    
    beth(i) = norm( Q(:, i+1) );            % Compute off-diagonal coefficients
    
    if ( abs(beth(i)) < sqrt(n)*eps ), break; end
    
    Q(:, i+1) = Q(:, i+1) / beth(i);        % Normalize
    
end

% i contains number of completed iterations

B = diag(aleph(1:i), 0) + diag(beth(1:(i-1)), +1) + diag(beth(1:(i-1)), -1);

[U, D] = cgal_eig(0.5*(B+B'));
[xi, ind] = min(diag(D));
v = Q(:, 1:i) * U(:, ind);

% Next lines are unnecessary in general, but I observed numerical errors in
% norm(v) at some experiments, so let's normalize it for robustness. 
nv = norm(v);
xi = xi*nv;
v = v/nv;
end


function [V,D] = cgal_eig(X)
% Eig in Lanczos based LMO solver sometimes fall into numerical issues. 
% This function replaces eig with a SVD based solver, in case eig does not
% converge. 
try
    [V,D]       = eig(X);
catch 
	warning('eig did not work. Using the svd based replacement instead.');
    [V,D,W]     = svd(X);
    d           = diag(D).' .* sign(real(dot(V,W,1)));
    [d,ind]     = sort(d);
    D           = diag(d);
    V           = V(:,ind);
end
end

%% ThinSVD Update
function [ U, S ] = updateThinSVDrank1( U1, S1, u2, s2 )
% FUNCTION:   [U, S] = updateThinSVD(...)
% PURPOSE:    Modification of given SVD under symmetric rank 1 perturbations.
%             [U, S, V] = svd(U1*S1*U1' + u2*s2*u2');
% REFERENCE:  Brand,M., Fast low-rank modification of thin singular value
%             decomposition, Linear Algebra and its Applications 415 (2006)
%             20-30, ELSEVIER, doi:10.1016/j.laa.2005.07.021

if ~any(u2), U = U1; S = S1; return; end

a = u2*sqrt(s2);
m = U1'*a;
P = a - U1*m;
Ra = norm(P);
if Ra, P = P/Ra; end

K = [S1, zeros(size(S1,1),1); zeros(1,size(S1,2)), 0] + [m; Ra]*[m; Ra]';

[u2, s2, ~] = svd(K,'econ');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Some heuristic rules to prune thin SVD to reduce memory consumption
while s2(end) <= max(diag(s2))*1e-9
    s2 = s2(1:end-1,1:end-1);
    u2(:,end) = [];
end
%               OR
% if s2(end) <= sum(diag(s2))*1e-9
%     s2 = s2(1:end-1,1:end-1);
%     u2(:,end) = [];
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

U = [U1, P]*u2;
S = s2;

end
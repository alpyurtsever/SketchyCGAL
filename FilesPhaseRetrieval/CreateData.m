%%  Create data for the Abstract Phase Retrieval experiments.
% Run this script to generate the random data for the experiment.

rng(0,'twister');

L = 10; % we get 10*n measurements

for n = [1e2,1e3,1e4,1e5,1e6]
    for MC = 1:20
        
        x_true = inf;                   % norm(x_true)^2 ~= 2*n with high probability
        while norm(x_true)^2 > 1.5*2*n  % but we add the sanity check to make sure that the problem is not ill-posed
            
            % Construct the true signal
            x_true = (randn(n,1) + 1i *randn(n,1));
            
            % Generate masks
            b1 = randsrc(n,L,[1i, -1i, 1, -1]);
            b2 = randsrc(n,L,[sqrt(3), 1/sqrt(2); 0.2, 0.8]);
            Masks = b1.*b2;
            clearvars b1 b2;
            
            % Take measurements
            b = evalAop_vec(x_true,Masks);
            
            % Save
            if ~exist(['./',num2str(L),'/data'],'dir'), mkdir(['./',num2str(L),'/data']); end
            save(['./',num2str(L),'/data/',num2str(n),'_',num2str(MC),'.mat']);
            
            fprintf('Data %d for n=%d generated... \n',MC,n);
            
        end
    end
end

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

%% Last edit: Alp Yurtsever - November 18, 2019
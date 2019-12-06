function out = spmult(i,j,s,n,x)
%MEXSPARSEMULT Computes "sparse(i,j,s,n,n)*x" without generating the sparse matrix

warning('You are running slow MATLAB implementation of SPMULT. Compile the mex codes for better performance.');

out = zeros(n,1);
for t = 1:length(s)
    IndI = i(t);
    IndJ = j(t);
    out(IndI) = out(IndI) + s(t)*x(IndJ);
end

end

%% Last edit: Alp Yurtsever - November 6, 2019

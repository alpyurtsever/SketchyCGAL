function out = mexSparseMult(i,j,s,n,x)
%MEXSPARSEMULT Computes "sparse(i,j,s,n,n)*x" without generating the sparse matrix

warning('You are running slow MATLAB implementation of MEXSPARSEMULT. Compile the mex codes for better performance.');

out = zeros(n,size(x,2));
for l = 1:size(x,2)
for t = 1:length(s)
    IndI = i(t);
    IndJ = j(t);
    out(l,IndI) = out(l,IndI) + s(t)*x(l,IndJ);
end
end

end

%% Last edit: Alp Yurtsever - January 19, 2020

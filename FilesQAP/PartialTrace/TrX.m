
function x = TrX(p,sys,dim)

% TRX   Partial trace
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    RHO = TrX(PSI,SYS,DIM) traces out the subsystems specified in
%    vector SYS of state PSI (a state vector or densitry matrix) whose
%    subsystem dimensions are specified by the vector DIM.


%% Copyright (C) 2004-2009 Toby Cubitt
%%
%% This program is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public License
%% as published by the Free Software Foundation; either version 2
%% of the License, or (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
%% MA 02110-1301, USA.


% check arguments
if any(sys > length(dim)) || any(sys < 0)
  error('Invalid subsystem in SYS')
end
if (length(dim) == 1 && mod(length(p)/dim,1) ~= 0)...
  || length(p) ~= prod(dim)
  error('Size of state PSI inconsistent with DIM');
end


% remove singleton dimensions
if exist('setdiff')
  % matlab
  sys = setdiff(sys,find(dim == 1));
else
  % octave
  sys = complement(find(dim == 1),sys);
end
dim = dim(find(dim ~= 1));


% calculate systems, dimensions, etc.
n = length(dim);
rdim = dim(end:-1:1);
keep = [1:n];
keep(sys) = [];
dimtrace = prod(dim(sys));
dimkeep = length(p)/dimtrace;


if any(size(p) == 1)
  % state vector
  if size(p,1) == 1
    p = p';
  end
  % reshape state vector to "reverse" ket on traced subsystems into a bra,
  % then take outer product
  perm = n+1-[keep(end:-1:1),sys];
  x = reshape(permute(reshape(p,rdim),perm),[dimkeep,dimtrace]);
  x = x*x';


else
  % density matrix

  % reshape density matrix into tensor with one row and one column index
  % for each subsystem, permute traced subsystem indices to the end,
  % reshape again so that first two indices are row and column
  % multi-indices for kept subsystems and third index is a flattened index
  % for traced subsystems, then sum third index over "diagonal" entries
  perm = n+1-[keep(end:-1:1),keep(end:-1:1)-n,sys,sys-n];
  x = reshape(permute(reshape(p,[rdim,rdim]),perm),...
              [dimkeep,dimkeep,dimtrace^2]);
  x = sum(x(:,:,[1:dimtrace+1:dimtrace^2]),3);

end
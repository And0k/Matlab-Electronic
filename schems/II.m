function z= II(varargin)
%all cell elements of varargin must be column vectors or scalars
z= cell2mat(varargin); n= numel(z);
N= int32(numel(varargin));
if n>N; z= reshape(z, n/N, N); end
z= 1./sum(1./z, 2);
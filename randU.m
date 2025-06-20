function [U]=randU(n,G)

% generates a random fuzzy partition of n subjects in g groups
%
%
U=rand(n,G);
U=repmat(1./sum(U,2),1,G).*U;
%U=bsxfun(@times, U, 1./(sum(U, 2)));
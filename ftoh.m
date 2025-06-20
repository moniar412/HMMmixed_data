function Uh=ftoh(Uf)
%
% converts a partition from fuzzy to hard
%
[n,nk]=size(Uf);
Uh=zeros(n,nk);
[a,ind]=max(Uf');
for i=1:n
    Uh(i,ind(i))=1;
end
function f=funct1(param,p0,cumrel1,O,th_idx,G,i)
%
%
f=0;
x=param(1);
mu=reshape(param(2:(1+O*G)),O,G)';
var=reshape(param((2+O*G):(1+O*G+O*G)),O,G)';

for g=1:G,
    f=f+p0(g)*normcdf(x,mu(g,th_idx(i)),sqrt(var(g,th_idx(i))));
end
f=(f-cumrel1(i))^2;
end



function f=funct(x,p0,mu,var,cumrel1,th_idx,G,i)
%
%
f=0;
for g=1:G,
    f=f+p0(g)*normcdf(x,mu(g,th_idx(i)),sqrt(var(g,th_idx(i))));
end
f=(f-cumrel1(i))^2;
end



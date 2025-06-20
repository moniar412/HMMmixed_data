function mri=mrand(N)
%
% modified rand index (Hubert & Arabie 1985, JCGS p.198)
%
n=sum(sum(N));
sumi=.5*(sum(sum(N').^2)-n);
sumj=.5*(sum(sum(N).^2)-n);
pb=2*sumi*sumj/(n*(n-1));
mri=(.5*(sum(sum(N.^2))-n)-pb)/((sumi+sumj)/2-pb);
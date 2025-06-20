function Gamm=Qm2Gamm(Qm)
%
% convert Qm parameters to Gamm
%
[G,G,Gm]=size(Qm);
Gamm=zeros(G,G,Gm);
one = ones(1,G);
for gm=1:Gm
    for g=1:G
        q = max(one*0.00001,min(Qm(g,:,gm),one*0.99999));
        Gamm(:,g,gm)=log(q/q(g));
    end
end
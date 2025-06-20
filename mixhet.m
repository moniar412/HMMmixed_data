function [U,Mmu,MOmeQ,MOmel,dif,like,awe,bic,aic,it]=mixhet(X,U,eps,dis)
% 
% heteroscedastic Gaussian mixture (ML solution)
%
% 13/06/2007
%
[n,J]=size(X);
G=size(U,2);
su=sum(U);
onesn=ones(n,1);
MOmeQ=zeros(J,G*J);
for g=1:G
    MOmeQ(:,(g-1)*J+(1:J))=eye(J);
end
MOmel=ones(J,G);
oJ=ones(J,1);
lfig=zeros(n,G);
p=(su/n)';
dif=1;
likeold=-Inf;
it=0;
%
while dif > eps,
   it=it+1;
   %
   % update Mmu
   Mmu = diag(1./su)*U'*X;
   %
   % update MOme and sig
   for g=1:G
           Sig=X-onesn*Mmu(g,:);
           Sig=Sig'*diag(U(:,g))*Sig/su(g);
       [P,L,Q]=svd(Sig);
       l=diag(L);
       MOmeQ(:,(g-1)*J+(1:J))=Q;
       MOmel(:,g)=min(oJ*10000,max(oJ*0.0001,l));
%        Rocci
%        MOmel(:,g)=min(oJ*10000,max(oJ*0.0001,l));
%        MOmel(:,g)=l;
   end
   %
   % update U
   for g=1:G
       Q=MOmeQ(:,(g-1)*J+(1:J));
       lam=MOmel(:,g);
       lfig(:,g)=-0.5*sum(log(lam))-0.5*sum(((X-onesn*Mmu(g,:))*Q*diag(1./sqrt(lam))).^2,2);
   end
   ind=(lfig<-7.0e2);
   U=exp(-7.0e2*ind+lfig.*(1-ind))*diag(p);
   U=repmat(1./sum(U,2),1,G).*U;
   su=sum(U);
   %
   % update p
   p=(su/n)';
   %
   % stopping rule
   warning off
   UlU=U.*log(U);
   ulu=UlU(:);
   ulu(find(isnan(ulu)))=[];
   sulu=sum(ulu);
   warning on
   like=n*sum(p.*log(p))+sum(sum(U.*lfig))-sulu;
   dif=(like-likeold);
   likeold=like;
   %
end
np=G-1+G*J+G*(J*J+J)*0.5;
%awe=2*like+2*sulu-2*(3/2+log(n))*np;
awe=1;
bic=2*like-log(n)*np;
aic=2*like-2*np;
if dis==1
%     disp(sprintf('mixhet(%g): dif=%g, iter=%g, like=%g, AWE=%g, BIC=%g, AIC=%g',G,dif,it,like,awe,bic,aic))
end
if dif<-0.000001
    disp('------------------------------------------- error -------------------------------------')
    dif
    MOmel
    10000*min(MOmel)
    sum(U)
end
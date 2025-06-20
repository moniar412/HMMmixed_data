function [f, plimix, plik]=plikObsAll(theta,m,th_idx,idxVar,X,Q,no_TH,nopairs,pairsindex)
%
% OUTPUT__________________________________
% the weighted likelihood of the whole mixture
% INPUT___________________________________
% param: a vector containing the parameters which need to be estimated.
% pg1: extra parameter, the component weights
% n_xixj: the vector of frequencies n_xixj_cicj
% G: the number of groups
% th_idx: the vector keeping track to which variable the thresholds belongs to
% T: the output of E.step
% P: the number of ordinal variables
%
% m & s: mean and variance continuous variables 
G=size(m,1);
P=size(X,2);
Xc=X(:,idxVar==0);
C=size(Xc,2);
n=size(Xc,1);
N=n;
Xo=X(:,idxVar==1);
O=size(Xo,2);



no.pairs=nopairs;
pairs.index=pairsindex;
    


V1=theta.V1;
V2=theta.V2;
choleg=theta.choleg;
cors=zeros(G,P*(P-1)/2);
sigma2=ones(G,P);
corsO=zeros(G,O*(O-1)/2);
% J=diag(chole.^2);
r=(find(tril(ones(size(ones(P,P))),-1)));
r0=(find(tril(ones(size(ones(O,O))),-1)));
idxO=find(idxVar==1);
SIGMA=ones(P,P,G);
for g=1:G,

    if g==1,
        S1=V1*V1';
        S2=V2*V2';
        Sig=S1+S2;
    else
    L=ones(Q,Q);
    L(tril(L,0)~=0)=choleg(g-1,:)';
    L=tril(L,0);
     S=L*L';
    Sig=V1*S*V1'+V2*V2';
    end
%     Sig=A*[S zeros(Q,P-Q); zeros(P-Q,Q) J]*A';
SIGMA(:,:,g)=Sig;
       sigma2(g,:)=diag(Sig)';
      c=tril(diag(sqrt(diag(Sig)).^-1)*Sig*diag(sqrt(diag(Sig)).^-1),-1);
    cors(g,:)=c(r)';
    Sig=Sig(idxO,idxO);

       c=tril(diag(sqrt(diag(Sig)).^-1)*Sig*diag(sqrt(diag(Sig)).^-1),-1);
    corsO(g,:)=c(r0)';

end


       




plimix=ones(N,G,O);

%   s2c=ones(G,O);
 theta.cors=corsO;

br1=threshUniv(theta,G,O,th_idx); 
t1=ones(N,size(Xo,2),G);
t2=ones(N,size(Xo,2),G);

for g=1:G,
    [gamma1_P,gamma1_Pext1,gamma1_Pext2]=vargamma(br1(g,:),th_idx,O);
    for p=1:size(Xo,2),
        val=Xo(:,p);
        % val=X(:,p);
        t1(:,p,g)=gamma1_Pext1{p}(val);
        t2(:,p,g)=gamma1_Pext2{p}(val);
    end
end 

idxC=find(idxVar==0);
idxO=find(idxVar==1);
for o=1:O,
   for g=1:G,
 s2c=SIGMA(idxO(o),idxO(o),g)-SIGMA(idxC,idxO(o),g)'*inv(SIGMA(idxC,idxC,g))*SIGMA(idxC,idxO(o),g);
% if s2c<0,
%     s2c=0.*05;
% end
for i=1:N, 
    muc=theta.mu(g,o)+SIGMA(idxC,idxO(o),g)'*inv(SIGMA(idxC,idxC,g))*(Xc(i,:)-m(g,:))';
    
    t3=(mvncdf(t2(i,o,g),t1(i,o,g),muc,s2c));
    t3(t3<eps)=eps;
%     t4=mvnpdf(Xc(i,mixed(2,o)),m(g,mixed(2,o)) ,SIGMA(mixed(2,o),mixed(2,o),g));
%      t4(t4<eps)=eps;
    plimix(i,g,o)=t3*mvnpdf(Xc(i,:),m(g,:),SIGMA(idxC,idxC,g));
end
      
    end
    

 
end

tables=Tables(Xo,O,no_TH);
n_xixj=frequency1(tables,O);
n=length(n_xixj);
TH=thresh(theta,th_idx,O,G);

% theta.cors(abs(theta.cors<eps))=0.05;
plik=ones(n,G);
for g=1:G,
    plik(:,g)=(PI_by_cell(TH(g,:),theta.cors(g,:),O,th_idx));
end







   
   ni=repmat(n_xixj',1,G);
% plik(plik<eps)=eps;
if O==2,
f=log(plik.*prod(plimix,3));
% sum(log(sum(U1,2)))+sum(sum(log(sum(plik,2))))+sum(sum(log(sum(plimix,2))));
% sum(sum(U.*U1-U.*log(U)))+sum(sum(T.*plik-T.*log(T)))+sum(sum(sum(Umix.*plimix-Umix.*log(Umix))));
 ind=(f<-7.0e2);
 f=exp(-7.0e2*ind+f.*(1-ind));

else
    
    
idP=cell(no.pairs,1);
idxP=ones(1,no.pairs);
for i=1:no.pairs,
idP{i}=[reshape(repmat(1:no_TH(pairs.index(1,i)),no_TH(pairs.index(2,i)),1)',1,[])'  reshape(repmat(1:no_TH(pairs.index(2,i)),no_TH(pairs.index(1,i)),1),1,[])'];
idxP(i)=size(idP{i},1);
end

CidxP=cumsum(idxP);
CidxP=[0 CidxP];
idPl=cell(no.pairs,1);
for i=1:no.pairs,
    idPl{i}=plik((CidxP(i)+1):CidxP(i+1),:);
end

Xoh=ones(N,G,no.pairs);
for i=1:no.pairs,
    for j=1:size(idP{i},1),
    idx=find(ismember(Xo(:,pairs.index(:,i)'),idP{i}(j,:),'rows'));
    Xoh(idx,:,i)=reshape(repmat(idPl{i}(j,:),length(idx),1)',G,[])';
    end
end
 plimix(plimix<eps)=eps;
 pr=prod((plimix),3);
  Xoh(Xoh<eps)=eps;
 pli=prod((Xoh),3);
        
    f=pr.*pli;
% 
 f=(f);

%   ind=(f<-7.0e2);
%   
%        f=exp(-7.0e2*ind+f.*(1-ind));

  
%        f=sum(log(sum(f,2)));

       
%        
%        
%    warning off      
% UlU=Umix.*log(Umix);    
%   ulu=UlU(:);
%    ulu(find(isnan(ulu)))=[];
%    sulu=sum(sum(ulu));     
%    warning on,
%    like1=c*N*sum(theta.pg.*log(theta.pg))+sum(sum(sum(Umix.*log(plimix))))-sulu;
%        
% 
%    
%     warning off      
% UlU=T.*log(T);    
%   ulu=UlU(:);
%    ulu(find(isnan(ulu)))=[];
%    sulu=sum(sum(ulu));     
%    warning on,
%    like2=length(n_xixj)*sum(theta.pg.*log(theta.pg))+sum(sum(T.*log(plik)))-sulu;
% 
%    
%    warning off      
% UlU=U.*log(U);    
%   ulu=UlU(:);
%    ulu(find(isnan(ulu)))=[];
%    sulu=sum(sum(ulu));     
%    warning on,
%    like3=N*sum(theta.pg.*log(theta.pg))+sum(sum(sum(U.*log(U1))))-sulu;
%        
% 
%    f=like1+like2+like3;
       
       %     sum(log(sum(U1,2)))+(sum(log(sum(plik,2)).*n_xixj'))+sum(sum(log(sum(plimix,2))));
%     sum(sum(U.*U1-U.*log(U)))+sum(sum(T.*plik.*ni-T.*log(T).*ni))+sum(sum(sum(Umix.*plimix-Umix.*log(Umix))));

%    
%     sum(log(sum(U1,2)))+sum(sum(log(sum(plik,2))))+sum(sum(log(sum(plimix,2))));
%     sum(log(sum((U1*diag(theta.pg)),2)))+sum(sum(n_xixj'.*log(sum((plik*diag(theta.pg)),2))))+sum(sum(log(sum((plimix*diag(theta.pg)),2))));

% f=sum(sum(U.*U1-U.*log(U)))+sum(sum(T.*plik.*ni-T.*log(T).*ni))+sum(sum(sum(Umix.*plimix-Umix.*log(Umix))));
end
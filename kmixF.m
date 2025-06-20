function [theta_init,Qm,Uplus]=kmixF(X,Y,no_TH,th_idx,P,G,Q,O,idxO,idxVar,pairsindex)
pairs.index=pairsindex;
no.pairs=O*(O-1)/2;
Yt=X;
U=randU(size(X,1),G);
T=size(Y,3);
Yc=Y(:,idxVar==0,:);
[U,Xmean,res]=km(X,U);
[out,Mmu,MOmeQ,MOmel,dif,lik,awe,bic,aic,it]=mixhet(X,U,10^-6,1);
p0=sum(out)/size(out,1);
% [B,ix] = sort(Mmu(:,1)');
% ix=ix';
% Mmu=Mmu(ix,:);
rst = randi([0 1000],1,5);
mu=cell(1,5);
mu1=cell(1,5);
SigA=cell(1,5);
Sig=cell(1,5);
Deltam=cell(1,5);
Deltam=cell(1,5);
Qm=cell(1,5);
Qm=cell(1,5);
U01=cell(1,5);
U0=cell(1,5);
Uplus1=cell(1,5);
Uplus=cell(1,5);
llk1=ones(1,5);
llk=ones(1,5);
for i=1:5,
[mu1{i},SigA{i},Deltam{i},Qm{i},U01{i},Uplus1{i},llk1(1,i),bic,iter]=mixhetgaus3(Y,G,10^-6,rst(i));
[mu{i},Sig{i},Deltam{i},Qm{i},U0{i},Uplus{i},llk(1,i),bic,iter]=mixhetgaus3(Yc,G,10^-6,rst(i));
end
id1=find(llk1==max(llk1));
id=find(llk==max(llk));
id=id(1);
mu=mu{id};
Qm=Qm{id};
out=Uplus{id};


p0=sum(Uplus{id})/size(Uplus{id},1);


[B,ix] = sort(mu(:,1),'descend');
mu=mu(ix,:);

Qm=Qm(ix,ix);

% var=var(ix,:);
p0=p0(ix);
out=out(:,ix);
Uplus=out;
outH=ftoh(out);


mu=Mmu;

mc=Mmu(:,idxVar==0);
[B,ix] = sort(mc(:,1),'descend');
mu=mu(ix,:);
% 
% Qm=Qm(ix,ix);
% 
% var=var(ix,:);
% p0=p0(ix);
% out=out(:,ix);
% Uplus=out;
% outH=ftoh(out);

variance=zeros(G,P);
R1=cell(G,1);

cors=zeros(G,no.pairs);
sigma=cell(G,1);
ch=cell(G,1);
v=cell(G,1);
r0=(find(tril(ones(size(ones(O,O))),-1)));

for g=1:G,
    sigma{g}=cov(X(find(outH(:,g)==1),:));
    R1{g}=corr(tiedrank(X(find(outH(:,g)==1),:)),'Type','Spearman'); 
    
    Ro=R1{g}(idxO,idxO);
    cors(g,:)=Ro(r0)';
    variance(g,:)=diag(sigma{g});
end

Xo=X(:,idxO);
c=unique(Xo);
h=hist(Xo,c);
var=variance;
%mixture thresholds
cumh=cumsum(h);
cumrel=bsxfun(@rdivide,cumh,cumh(size(cumh,1),:));
cumrel(size(cumrel,1),:)=[];
cumrel1=reshape(cumrel,1,[]);
th=cell(1,O);
for j=1:O,
    c=unique(Xo(:,j));
    h=hist(Xo(:,j),c)';
    
    %mixture thresholds
    cumh=cumsum(h);
    cumrel=bsxfun(@rdivide,cumh,cumh(size(cumh,1),:));
    if length(cumrel)==length(th_idx(th_idx==j)),
        cumrel(size(cumrel,1),:)=[];
        cumrel=[cumrel; cumrel(size(cumrel,1))+rand(1)];
        cumrel1=reshape(cumrel,1,[]);
    elseif length(cumrel)<length(th_idx(th_idx==j)),
        cumrel(size(cumrel,1),:)=[];
        cumrel=[cumrel cumrel(size(cumrel,1))+rand(length(th_idx(th_idx==j))-size(cumrel,1))];
        cumrel1=reshape(cumrel,1,[]);
    else
        cumrel(size(cumrel,1),:)=[];
        cumrel1=reshape(cumrel,1,[]);
    end
    
    % fun=cell(length(th_idx),1);
    th{j}=zeros(1,length(th_idx(th_idx==j)));
    x0=normrnd(0,1,1,length(th_idx(th_idx==j)));
    % x0=randi([-10 10]);
    for i=1:length(th_idx(th_idx==j)),
                th{j}(i)=fminsearch(@(x)funct(x,p0,mu(:,idxO(j)),var(:,idxO(j)),cumrel1,th_idx,G,i),x0(i));
    end
  th{j}=  sort(th{j});
 end
th=cell2mat(th);

gamma=th;
for g=1:G,
    for i=1:O,
    ga=th(th_idx==i);   
    ga=(ga-ga(1));
    ga=ga/ga(2);
    gamma(th_idx==i)=ga;
      end
end

theta=struct('pg',p0,'mu',mu(:,idxO),'gamma',gamma,'sigma',ones(1,O),'cors',cors);
param=mu(:,idxO);
param=param(:)';
%interior-point
%sqp
%active-set
options=optimset('fmincon');
options = optimset(options,'Algorithm','interior-point','Display','on','LargeScale','on','TolFun',1.0e-6,'TolCon',1.0e-6,'MaxFunEvals',50*length(param),'MaxIter',10*length(param),'UseParallel','Always');
lb=repmat(-10,1,length(param));
ub=repmat(+10,1,length(param));


Uplusplus=repmat(out,1,1,size(pairs.index,2));



idP=cell(1,size(pairs.index,2));
dimPair=ones(1,size(pairs.index,2));
for idx=1:size(pairs.index,2),
    id1=repmat(1:no_TH(pairs.index(1,idx)),1,no_TH(pairs.index(2,idx)),1)';
    id2=repmat(1:no_TH(pairs.index(2,idx)),no_TH(pairs.index(1,idx)),1);
    id2=id2(:);
    idP{idx}=[id1,id2];
    dimPair(idx)=size(idP{idx},1);
end

cL=cumsum(dimPair);
cLU=[cL(1:(size(pairs.index,2)-1)) cL((size(pairs.index,2)))];
cLL=[1 cL(1:(size(pairs.index,2)-1))+1];
ct=(1:T)*size(Y,1);
ctl=[1 ct(1:(T-1))+1];
ctu=[ct(1:(T-1)) ct(T)];

% Uplust=Uplusplus(1:size(Y,1),:,:);
% X=Y(:,:,1);
% tic,[par,f]=fmincon(@(param)funct1opt(param,theta,X,O,th_idx,idxVar,G,pairs.index,Uplusplus),param,[],[],[],[],lb,ub,[],options);
% toc
% tic,
[par,f]=fmincon(@(param)funct1opt2(param,theta,Y,O,th_idx,G,idxVar,pairs.index,Uplusplus,idP,cLU,cLL,ctl,ctu,dimPair),param,[],[],[],[],lb,ub,[],options);
% toc
% tic,[par,f]=fmincon(@(param)funct1optMEX(param,theta,X,O,th_idx,idxVar,G,pairs.index,Uplusplus),param,[],[],[],[],lb,ub,[],options);
% toc
muh=reshape(par,G,O);
[B,ix] = sort(muh(:,1),'ascend');
ix=ix';
mu(:,idxO)=muh(ix,:);
theta.mu=muh;


for o=1:O,
ido=th(th_idx==o);
var(1,idxO(o))=(1/(ido(2)-ido(1)))^2;
end

param=var(:,idxO);
param=param(:)';
lb=repmat(0.05,1,length(param));
ub=repmat(10,1,length(param));
options = optimset(options,'Algorithm','interior-point','Display','on','LargeScale','on','TolFun',1.0e-6,'TolCon',1.0e-6,'MaxFunEvals',50*length(param),'MaxIter',10*length(param),'UseParallel','Always');

% tic,
% [par,f]=fmincon(@(param)funct1opt2s(param,theta,Y,O,th_idx,G,idxVar,pairs.index,Uplusplus,idP,cLU,cLL,ctl,ctu,dimPair),param,[],[],[],[],lb,ub,[],options);
%  toc
%sigma2=reshape(par,G,O);
 %var(:,idxO)=sigma2;

 %for g=1:G,
%[c1 s1 l1]=pca(X(find(outH(:,g)==1),:));
% fg=diag(c1*diag(l1)*c1');
% var(g,idxO)=fg(idxO);
%end
%var(:,idxO)=sigma2;
Sigh=ones(P,P,G);
for g=1:G,
Sigh(:,:,g)=diag(sqrt((var(g,:))))*R1{g}*diag(sqrt((var(g,:))));
end



r0=(find(tril(ones(size(ones(O,O))),-1)));
r=(find(tril(ones(size(ones(P,P))),-1)));
ch=zeros(G,P*(P+1)/2);
for g=1:G,
    S=Sigh(:,:,g)';
     [e d]=eig(S);
     d(d<0)=0.01;
     S=e*d*e';
      v=chol(S);
      v=v';
    ch(g,:)=v(find(tril(ones(size(v)),0)))';
    %V*[S zeros(Q,P-Q); zeros(P-Q,Q) eye(P-Q)]*V';
   
end
    
P=size(X,2);
Xc=X(:,idxVar==0);
C=size(Xc,2);
N=size(Xc,1);
oP=ones(P,1);
Xo=X(:,idxVar==1);
O=size(Xo,2);
no.pairs=O*(O-1)/2;
j=find(gamma>1);
m=(reshape(mu(:,idxVar==1)',1,G*O));
muc=(reshape(mu(:,idxVar==0)',1,G*(C)));

param=[reshape(ch',1,[]) ] ;

costL=cell(P,1);costU=cell(P,1);
for i=1:P
     %costL{i}=[0.05,repmat(-5,1,P-i)];
    costL{i}=[-5,repmat(-5,1,P-i)];
    costU{i}=[5,repmat(5,1,P-i)];
end
options=optimset('fmincon');
options = optimset(options,'Algorithm','interior-point','Display','on','LargeScale','on','TolFun',1.0e-6,'TolCon',1.0e-6,'MaxFunEvals',30*length(param),'MaxIter',5*length(param),'UseParallel','Always');
lb=[  repmat(cell2mat(costL'),1,G) ];
ub=[ repmat(cell2mat(costU'),1,G) ];
%
%param,br,mu,m,th_idx,P,O,G,idxVar,Uplusplus,pairsindex,Yt
%
options = optimset(options,'Algorithm','sqp','Display','on','LargeScale','on','TolFun',1.0e-6,'TolCon',1.0e-6,'MaxFunEvals',500*length(param),'MaxIter',50*length(param),'UseParallel','Always');

[par,f]=fmincon(@(param)myfunmixN0(param,gamma,mu,mu(:,idxVar==0),th_idx,P,O,G,idxVar,Uplusplus,pairsindex,Yt),param,[],[],[],[],lb,ub,[],options);

a2=reshape(par,P*(P+1)/2,G)';

for g=1:G
    L=ones(P,P);
    L(tril(L,0)~=0)=a2(g,:)';
    L=tril(L,0);
    Sigh(:,:,g)=L*L';
end

%Sigh=Sigh(:,:,ix);
A=chol((Sigh(:,:,1)))';
%fare pca sulla matrice between

[coeff,score,latent] = pca(Yt);

%[coeff,score,latent] = pca(Yt(outH(:,1)==1,:));

[rr,ss]=sort(latent,'descend');
A=A(:,ss);
B=coeff*diag(sqrt(latent));
A1=A(:,1:Q);
A2=A(:,Q+1:P);
[Qs,Rs]=qr(A1(1:Q,1:Q)');
V1=A1*Qs;
[Qs2,Rs2]=qr(A2(1:(P-Q),:)');
V2=A2*Qs2;

V=[V1 V2];
%A1=A(:,1:Q);
%
%non è generalizzabile
%A2=[A(((P-Q)+1):P,((P-Q)+1):P);A(((P-Q)+1):P,1:(P-Q))];
%[QQ,RR]=qr(A1(1:Q,:)');
%V1=A1*(1*QQ);
%[QN,RN]=qr(A2(1:(P-Q),:)');
%V2=(A2*QN);

% A1=A(:,1:Q);
% A2=A(:,Q+1:P);
% [Qs,Rs]=qr(A1(1:Q,1:Q)');
% V1=A1*Qs;
% 
% [Qs2,Rs2]=qr(A2(Q+1,:)');
% V2=A2*Qs2;
% 
% V=[V1 V2];
% V=[V1 V2];
etah=ones(G,P);
for g=1:G,
   etah(g,:)=(pinv(V)*mu(g,:)')';

end
etah(:,ss((Q+1):P))=repmat(mean(etah(:,ss((Q+1):P)),1),G,1);



%Sigh1=Sigh;
r0=(find(tril(ones(size(ones(O,O))),-1)));
r=(find(tril(ones(size(ones(P,P))),-1)));
ch=zeros(G,Q*(Q+1)/2);
for g=1:G,
    S=inv(A)*(Sigh(:,:,g))*inv(A)';
    S=S(1:Q,1:Q);
     [e d]=eig(S);
     %d(d<0)=0.1;
     %S=e*d*e';
     S=S+eye(Q)+0.01;
      v=chol(S);
      v=v';
    ch(g,:)=v(find(tril(ones(size(v)),0)))';
    %V*[S zeros(Q,P-Q); zeros(P-Q,Q) eye(P-Q)]*V';
   
end
    




theta_init=struct('mu',mu,'pg',p0,'eta',etah,'choleg',ch(2:end,:),'cors',cors,'V1',V1,'V2',V2,'gamma',gamma,'sigma',var(:,idxO));
end
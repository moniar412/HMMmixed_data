function [theta, MUc,f]=M_stepMixed(Uplusplus,theta,Deltam,nopairs,pairsindex,G,Q,Yt,idxVar,th_idx)

% OUTPUT__________________________________
% compute the M-step; the main output is theta, i.e. a list of the updated parameters - pg, means,variances,
% correlations and gamma's. Since there is not a closed form for the parameter estimates, an iterative optimization
% is used (fminsearch). By default fminsearch performs minimization, but it maximize because the function is
% multiplied by -1.
% INPUT___________________________________
% theta: a list containing all the parameters - pg, means, variances, correlations,
% gamma's (the invariant thresholds) by group.
% n_xixj: the vector of frequencies n_xixj_cicj
% G: the number of groups
% th_idx: a vector keeping track to which variable the thresholds belongs to
% P: the number of ordinal variables
% no.pairs: the number of pairs of variables.
%
X=Yt;
P=size(X,2);
Xc=X(:,idxVar==0);
N=size(Xc,1);
Xo=X(:,idxVar==1);
O=size(Xo,2);
no.pairs=O*(O-1)/2;
% t=length(th_idx);
% to obtain the inputs for the log-likelihood
j=find(theta.gamma>1);
ColN=cell(1,P-Q);
Col=cell(1,Q);
oP=ones(P,1);
for q=1:Q
    Col{q}=repmat(q,1,P-q+1);
end
ColIdx=cell2mat(Col);
for p=1:P-Q,
    ColN{p}=repmat(p,1,P-Q-p+1);
end
ColIdxN=cell2mat(ColN);
V1=theta.V1;
V2=theta.V2;
for i=1:Q
    idx=1:length(ColIdx(ColIdx==i));
    if length (idx)<P,
        idx=idx+ (P-length(idx));
    else idx=idx;
    end
    t0=V1(idx,i)';
    if i==1,
        Vid=t0;
    else
        Vid=[Vid t0];
    end
end
for i=1:P-Q
    idx=1:length(ColIdxN(ColIdxN==i));
    if length (idx)<P-Q,
        idx=idx+ (P-Q-length(idx));
    else idx=idx;
    end
    t0=V2(idx,i)';
    if i==1,
        Vid2=t0;
    else
        Vid2=[Vid2 t0];
    end
end
V21=V2(P-Q+1:end,:);
V21=V21(:)';
V2=V2(1:Q,:);
VidT=(Vid);
Vid2T=(Vid2);
V21T=(V21);
O=length(find(idxVar==1));
constr=cell(1,Q);
for i=1:Q,
    constr{i}=[-3 repmat(-3,1,(Q-i))];
   % constr{i}=[0.05 repmat(e0,1,(Q-i))];
end
constr=cell2mat(constr);
param=[reshape(theta.eta(:,1:Q)',1,[]),VidT,reshape(theta.choleg',1,[]),V21T,Vid2T ,reshape(theta.eta(1,Q+1:end)',1,[])];
lb=zeros(1,length(param)); ub=zeros(1,length(param));
lb=[repmat(-10,1,Q*G) repmat(-5,1,(Q*(P-Q)+Q*(Q+1)/2)) repmat(constr,1,G-1) repmat(-5,1,(Q*(P-Q)+(P-Q)*(P-Q+1)/2))   repmat(-10,1,P-Q)];
ub=[repmat(10,1,Q*G) repmat(5,1,(Q*(P-Q)+Q*(Q+1)/2)) repmat(5,1,(G-1)*Q*(Q+1)/2) repmat(5,1,(Q*(P-Q)+(P-Q)*(P-Q+1)/2))  repmat(10,1,P-Q)];


param=[param theta.alpha];
ub=[ub repmat(10,1,length(theta.alpha))];
lb=[lb repmat(0,1,length(theta.alpha))];


options=optimset('fmincon');

%PER RENDERLO PI  VELOCE
options = optimset(options,'Algorithm','interior-point','Display','on','LargeScale','on','UseParallel','Always','TolFun',1.0e-6,'TolCon',1.0e-6,'MaxFunEvals',10*(length(param)),'MaxIter',10);


[par,f]=fmincon(@(param)myfunmixedC(param,th_idx,P,G,Q,j,ColIdx,ColIdxN,idxVar,Uplusplus,nopairs,pairsindex,N,Yt,Deltam'),param,[],[],[],[],lb,ub,[],options);

 param=par;
O=length(find(idxVar==1));
eta1=(reshape(param(1:Q*G),Q,G))';
V1c=param(Q*G+1:Q*G+(Q*(Q+1)/2+Q*(P-Q)));
ch=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2);
V21=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q));
V2c=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+((P-Q)*(P-Q+1)/2));
eta2=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+((P-Q)*(P-Q+1)/2)+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+((P-Q)*(P-Q+1)/2)+P-Q);
b=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+((P-Q)*(P-Q+1)/2)+P-Q+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+((P-Q)*(P-Q+1)/2)+P-Q+length(th_idx)-2*O);
a=1:length(th_idx);
a(j)=[];  bb=repmat([0 1],1,O);
br=zeros(1,length(th_idx)); br(j)=b; br(a)=bb;
eta=[eta1';repmat(eta2,G,1)']';
ColIdxT=repmat(ColIdx,1);
V1=zeros(P,Q);
V2=zeros(P,P-Q);
for q=1:Q,
    v1=vec2mat( V1c(ColIdxT==q),P-(q-1))';
    V1(:,q)=[zeros(1,q-1) v1']';
end
V22=reshape(V21,Q,P-Q);
V24=zeros(P-Q,P-Q);
ColIdxNT=repmat(ColIdxN,1);
for p=1:P-Q,
    v1=vec2mat(V2c(ColIdxNT==p),P-Q-(p-1))';
    V24(:,p)=[zeros(1,p-1) v1']';
end
V2=[V24;V22];
choleg=reshape(ch,Q*(Q+1)/2,G-1)';
V=[V1 V2];
%V=eye(P);

mu=ones(G,P);
for g=1:G,
    mu(g,:)=(V*eta(g,:)')';
end
% eta1=(V1*eta1')';
% mu=[eta1 ];
% eps=10^-2;
alpha=b;
for i=1:O,
     ga=cumsum(br(th_idx==i));
      
     br(th_idx==i)=ga;
end
cors=zeros(G,P*(P-1)/2);
sigma2=ones(G,P);
corsO=zeros(G,O*(O-1)/2);
corsP=zeros(G,P*(P-1)/2);
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
    %[P0,L,Q0]=svd(Sig);
    %l=diag(L);
    %l=min(oP*10000,max(oP*0.0001,l));
    %Sig=Q0*diag(l)*Q0';
    
    SIGMA(:,:,g)=Sig;
    sigma2(g,:)=diag(Sig)';
    c=tril(diag(sqrt(diag(Sig)).^-1)*Sig*diag(sqrt(diag(Sig)).^-1),-1);
    cors(g,:)=c(r)';
    Sig=Sig(idxO,idxO);
    c=tril(diag(sqrt(diag(Sig)).^-1)*Sig*diag(sqrt(diag(Sig)).^-1),-1);
    corsO(g,:)=c(r0)';
end
m=mu(:,idxVar==0,:);
theta1=struct('pg',Deltam','mu',mu,'eta',eta,'sigma',sigma2,'alpha',alpha,'V1',V1,'V2',V2,'choleg',choleg,'cors',cors,'corsO',corsO,'gamma',br,'SIGMA',SIGMA);
theta1.mu=mu(:,idxVar==1);
theta1.sigma=sigma2(:,idxVar==1);
theta1.cors=corsO;
theta=theta1;
m0=m;
MUc=m0;
end



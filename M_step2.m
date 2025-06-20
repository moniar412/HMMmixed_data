function [theta, MUc,f]=M_step2(T,gamma0,gamma0p,n_xixj,theta,delta,nopairs,pairsindex,no_TH,G,Q,Y,idxVar,th_idx,n)

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
X=Y;
P=size(X,2);
Xc=X(:,idxVar==0);
C=size(Xc,2);
N=size(Xc,1);

Xo=X(:,idxVar==1);
O=size(Xo,2);

no.pairs=O*(O-1)/2;
% t=length(th_idx);
% to obtain the inputs for the log-likelihood
j=find(theta{1}.gamma>1);

ColN=cell(1,P-Q);
Col=cell(1,Q);

for q=1:Q
Col{q}=repmat(q,1,P-q+1);
end
ColIdx=cell2mat(Col);
for p=1:P-Q,
    ColN{p}=repmat(p,1,P-Q-p+1);
end
ColIdxN=cell2mat(ColN);







Vid=cell(1,T);
V21=cell(1,T);
Vid2=cell(1,T);

for t=1:T,
    V1=theta{t}.V1;
V2=theta{t}.V2;
for i=1:Q
       idx=1:length(ColIdx(ColIdx==i));
    if length (idx)<P,
        idx=idx+ (P-length(idx));
    else idx=idx;
    end
    t0=V1(idx,i)';
    if i==1,
        Vid{t}=t0;
    else
    Vid{t}=[Vid{t} t0];
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
        Vid2{t}=t0;
    else
    Vid2{t}=[Vid2{t} t0];
    end
end

V21{t}=V2(P-Q+1:end,:);
V21{t}=V21{t}(:)';
V2=V2(1:Q,:);
    
end

VidT=cell2mat(Vid);
Vid2T=cell2mat(Vid2);
V21T=cell2mat(V21);

O=length(find(idxVar==1));

constr=cell(1,Q);
 for i=1:Q,
     constr{i}=[0.05 repmat(-5,1,(Q-i))];
 end
constr=cell2mat(constr);
param=[reshape(theta{1}.eta(:,1:Q)',1,[]),VidT,reshape(theta{1}.choleg',1,[]),V21T,Vid2T ,reshape(theta{1}.eta(1,Q+1:end)',1,[])];

lb=zeros(1,length(param)); ub=zeros(1,length(param));
lb=[repmat(-20,1,G*Q) repmat(-10,1,(Q*(P-Q)+Q*(Q+1)/2)*T) repmat(constr,1,G-1) repmat(-10,1,(Q*(P-Q)+(P-Q)*(P-Q+1)/2)*T)   repmat(-20,1,P-Q)];
%  lb=[repmat(-Inf,1,G*P+P*P) repmat(-5,1,G*Q*(Q+1)/2) repmat(-5,1,P-Q)  repmat(1,1,length(th_idx)-2*P)];
ub=[repmat(+20,1,G*Q) repmat(10,1,(Q*(P-Q)+Q*(Q+1)/2)*T) repmat(5,1,(G-1)*Q*(Q+1)/2) repmat(10,1,(Q*(P-Q)+(P-Q)*(P-Q+1)/2)*T)  repmat(+20,1,P-Q)];


param=[param theta{1}.gamma(theta{1}.gamma>1)];
ub=[ub repmat(+20,1,length(th_idx)-2*O)];
lb=[lb repmat(+1,1,length(th_idx)-2*O)];
options=optimset('fmincon');
options = optimset(options,'Algorithm','interior-point','LargeScale','off', 'TolFun',1.0e-6,'TolCon',1.0e-6,'MaxFunEvals',100,'MaxIter',10,'Display','off');
[par,f]=fmincon(@(param)myfun4(param,th_idx,nopairs,pairsindex,no_TH,P,G,Q,j,n,ColIdx,ColIdxN,idxVar,gamma0,gamma0p,n_xixj,N,Y,T,delta),param,[],[],[],[],lb,ub,[],options);


param=par;
O=length(find(idxVar==1));
C=P-O;
Xc=X(:,idxVar==0);

Xo=X(:,idxVar==1);



eta1=(reshape(param(1:Q*G),Q,G))';
V1c=param(Q*G+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))*T);
ch=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))*T+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))*T+(G-1)*Q*(Q+1)/2);


V21=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))*T+(G-1)*Q*(Q+1)/2+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))*T+(G-1)*Q*(Q+1)/2+Q*(P-Q)*T);
V2c=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))*T+(G-1)*Q*(Q+1)/2+Q*(P-Q)*T+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))*T+(G-1)*Q*(Q+1)/2+Q*(P-Q)*T+((P-Q)*(P-Q+1)/2)*T);


eta2=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))*T+(G-1)*Q*(Q+1)/2+Q*(P-Q)*T+((P-Q)*(P-Q+1)/2)*T+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))*T+(G-1)*Q*(Q+1)/2+Q*(P-Q)*T+((P-Q)*(P-Q+1)/2)*T+P-Q);

b=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))*T+(G-1)*Q*(Q+1)/2+Q*(P-Q)*T+((P-Q)*(P-Q+1)/2)*T+P-Q+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))*T+(G-1)*Q*(Q+1)/2+Q*(P-Q)*T+((P-Q)*(P-Q+1)/2)*T+P-Q+length(th_idx)-2*O);






a=1:length(th_idx);
a(j)=[];  bb=repmat([0 1],1,O);
br=zeros(1,length(th_idx)); br(j)=b; br(a)=bb;

eta=[eta1';repmat(eta2,G,1)']';
ColIdxT=repmat(ColIdx,1,T);
V1=zeros(P,Q,T);
V2=zeros(P,P-Q,T);


for q=1:Q,
v1=vec2mat( V1c(ColIdxT==q),P-(q-1))';
for t=1:T,
V1(:,q,t)=[zeros(1,q-1) v1(:,t)']';
end
end

V22=reshape(V21,Q,P-Q,T);
V24=zeros(P-Q,P-Q,T);




ColIdxNT=repmat(ColIdxN,1,T);
for p=1:P-Q,
    
    v1=vec2mat(V2c(ColIdxNT==p),P-Q-(p-1))';
for t=1:T,
    V24(:,p,t)=[zeros(1,p-1) v1(:,t)']';

end
end

V2=[V24;V22];

choleg=reshape(ch,Q*(Q+1)/2,G-1)';

V=[V1 V2];

mu=ones(G,P,T);
for t=1:T,
for g=1:G,
    mu(g,:,t)=(V(:,:,t)*eta(g,:)')';
end
end
% eta1=(V1*eta1')';

% mu=[eta1 ];
% eps=10^-2;
for i=1:O,
    
    br(th_idx==i)=sort(br(th_idx==i));
   ga=br(th_idx==i);
%    if ga(end)-ga(end-1)<eps,
%        ga(end)=ga(end)+2*eps;
%    end
   br(th_idx==i)=ga;
end

cors=zeros(G,P*(P-1)/2,T);
sigma2=ones(G,P,T);
corsO=zeros(G,O*(O-1)/2,T);
corsP=zeros(G,P*(P-1)/2,T);
% J=diag(chole.^2);
r=(find(tril(ones(size(ones(P,P))),-1)));
r0=(find(tril(ones(size(ones(O,O))),-1)));
idxO=find(idxVar==1);
SIGMA=ones(P,P,G,T);

for t=1:T,
for g=1:G,

    if g==1,
        S1=V1(:,:,t)*V1(:,:,t)';
        S2=V2(:,:,t)*V2(:,:,t)';
        Sig=S1+S2;
    else
    L=ones(Q,Q);
    L(tril(L,0)~=0)=choleg(g-1,:)';
    L=tril(L,0);
     S=L*L';
    Sig=V1(:,:,t)*S*V1(:,:,t)'+V2(:,:,t)*V2(:,:,t)';
    end
%     Sig=A*[S zeros(Q,P-Q); zeros(P-Q,Q) J]*A';
SIGMA(:,:,g,t)=Sig;
       sigma2(g,:,t)=diag(Sig)';
      c=tril(diag(sqrt(diag(Sig)).^-1)*Sig*diag(sqrt(diag(Sig)).^-1),-1);
    cors(g,:,t)=c(r)';
    Sig=Sig(idxO,idxO);

       c=tril(diag(sqrt(diag(Sig)).^-1)*Sig*diag(sqrt(diag(Sig)).^-1),-1);
    corsO(g,:,t)=c(r0)';

end
end
m=mu(:,idxVar==0,:);



theta=cell(1,T);
m0=cell(1,T);
for t=1:T,
theta1=struct('pg',delta,'mu',mu(:,:,t),'eta',eta,'sigma',sigma2(:,:,t),'V1',V1(:,:,t),'V2',V2(:,:,t),'choleg',choleg,'cors',cors(:,:,t),'corsO',corsO(:,:,t),'gamma',br,'SIGMA',SIGMA(:,:,:,t));
theta1.mu=mu(:,idxVar==1);
theta1.sigma=sigma2(:,idxVar==1);
theta1.cors=corsO;

theta{t}=theta1;
m0{t}=m(:,:,t);
end




MUc=m0;


end



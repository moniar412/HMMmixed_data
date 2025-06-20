function f= myfunmixedC(param,th_idx,P,G,Q,j,ColIdx,ColIdxN,idxVar,Uplusplus,nopairs,pairsindex,N,Yt,delta)
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
%
%#codegen
pairs.index=pairsindex;
O=length(find(idxVar==1));
Xc=Yt(:,idxVar==0);
Xo=Yt(:,idxVar==1);
oP=ones(P,1);
eta1=(reshape(param(1:Q*G),Q,G))';
V1c=param(Q*G+1:Q*G+(Q*(Q+1)/2+Q*(P-Q)));
ch=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2);
V21=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q));
V2c=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+((P-Q)*(P-Q+1)/2));
eta2=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+((P-Q)*(P-Q+1)/2)+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+((P-Q)*(P-Q+1)/2)+P-Q);
gamma=param(Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+((P-Q)*(P-Q+1)/2)+P-Q+1:Q*G+(Q*(Q+1)/2+Q*(P-Q))+(G-1)*Q*(Q+1)/2+Q*(P-Q)+((P-Q)*(P-Q+1)/2)+P-Q+length(th_idx)-2*O);
eta=[eta1';repmat(eta2,G,1)']';
ColIdxT=repmat(ColIdx,1);
V1=zeros(P,Q);

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

V1=V(:,1:Q);
V2=V(:,Q+1:P);
mu=ones(G,P);
for g=1:G,
    mu(g,:)=(V*eta(g,:)')';
end
% eta1=(V1*eta1')';
% mu=[eta1 ];
% eps=10^-2;

a=1:length(th_idx);
idb=ismember(1:length(th_idx),j);
if isempty(j)==1
    br=[0 1 0 1];
else
    a(j)=[];
    bb=repmat([0 1],1,O);
    br=zeros(1,length(th_idx));
    br(j)=gamma;
    br(find(idb==0))=bb;
    for i=1:Q
        ga=cumsum(br(th_idx==i));
        br(th_idx==i)=ga;
    end
end

%cors=zeros(G,P*(P-1)/2);
sigma2=ones(G,P);
corsO=zeros(G,O*(O-1)/2);
% J=diag(chole.^2);
r=(find(tril(ones(size(ones(P,P))),-1)));
r0=(find(tril(ones(size(ones(O,O))),-1)));
idxO=find(idxVar==1);
SIGMA=ones(P,P,G);
for g=1:G,
    if g==1,
        Sig=V1*V1'+V2*V2';
    else
        L=ones(Q);
        L(tril(L)~=0)=choleg(g-1,:)';
        L=tril(L);
        Sig=V1*(L*L')*V1'+V2*V2';
    end
    
    
    SIGMA(:,:,g)=Sig;
end
m=mu(:,idxVar==0,:);
m0=m;
plimix=ones(N,G,size(pairs.index,2));
%   s2c=ones(G,O);
br1=repmat(br,G,1);
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
t11=zeros(N,P,G);
t22=zeros(N,P,G);
idxC=find(idxVar==0);
idxO=find(idxVar==1);
t11(:,idxO,:)=t1;
t22(:,idxO,:)=t2;
mu0=zeros(G,P);
mu0(:,idxO)=mu(:,idxVar==1,:);
parfor idx=1:size(pairs.index,2),
    for g=1:G,
        s2c1=SIGMA(idxC,idxC,g);
        S2=s2c1+0.001*eye(length(idxC));
        %iS2=inv(S2);
        Rob=SIGMA(idxO(pairs.index(:,idx)),idxC,g)/S2;
        s2c=SIGMA(idxO(pairs.index(:,idx)),idxO(pairs.index(:,idx)),g)-Rob*SIGMA(idxC,idxO(pairs.index(:,idx)),g);
        %s2c
        s2c=s2c+0.0001*eye(2);
        %[P0,L,Q0]=svd(s2c);
        %l=diag(L);
        %l=min(10000,max(0.0001,l));
        %s2c=Q0*diag(l)*Q0';
        [s2,r2c]=cov2corr(s2c);
        muc=mu0(g,idxO(pairs.index(:,idx)))'+SIGMA(idxC,idxO(pairs.index(:,idx)),g)'/S2*(Xc-m(g,:))';
        low=(t11(:,idxO(pairs.index(:,idx)),g)-muc')./repmat(sqrt(diag(s2c)'),size(t11,1),1);
        up=(t22(:,idxO(pairs.index(:,idx)),g)-muc')./repmat(sqrt(diag(s2c)'),size(t11,1),1);
        B=[low,up];
        %B=round(B,5);
        [uu,idu,mm]=unique( B, 'rows');
        t3=bvnP(B(:,1),B(:,3),B(:,2),B(:,4),r2c(1,2));
        t4=mvnpdf(Xc,m(g,:),S2);
        %t4(t4<eps)=eps;
        plimix(:,g,idx)=t3(mm).*t4;
       
    end
end
f=log(plimix);
%ind=(f<log(eps));
%f=(log(eps)*ind+f.*(1-ind));
f=f.*Uplusplus;
f=-sum(f(:));
end

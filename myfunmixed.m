function f= myfunmixed(param,th_idx,P,G,Q,j,ColIdx,ColIdxN,idxVar,Uplusplus,nopairs,pairsindex,N,Yt,delta)
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
mu=ones(G,P);
for g=1:G,
    mu(g,:)=(V*eta(g,:)')';
end
% eta1=(V1*eta1')';
% mu=[eta1 ];
% eps=10^-2;

a=1:length(th_idx);
if isempty(j)==1
    br=[0 1 0 1];
else
    a(j)=[]; 
    bb=repmat([0 1],1,Q);
    br=zeros(1,length(th_idx));
    br(j)=gamma;
    br(a)=bb;
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
      
       %[Pc,Lc,Qc]=svd(Sig);
%         l=diag(Lc);
%         l=min(10000,max(0.0001,l));
%         Sig=Qc*diag(l)*Qc';
   Sig=Sig+ 0.0001*eye(P); 
       SIGMA(:,:,g)=Sig;
end
m=mu(:,idxVar==0,:);
m0=m;
plimix=ones(N,G,size(pairs.index,2));
%   s2c=ones(G,O);
br1=threshUniv(theta1,G,O,th_idx);
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
mu0(:,idxO)=theta1.mu;
for idx=1:size(pairs.index,2),
    for g=1:G,
        Rob=SIGMA(idxO(pairs.index(:,idx)),idxC,g)/SIGMA(idxC,idxC,g);
        s2c=SIGMA(idxO(pairs.index(:,idx)),idxO(pairs.index(:,idx)),g)-Rob*SIGMA(idxC,idxO(pairs.index(:,idx)),g);
        s2c=s2c+0.01*eye(2);
      
        
        
        % if s2c<0,
        %     s2c=0.*05;
        % end
        for n=1:N,
            muc=mu0(g,idxO(pairs.index(:,idx)))'+Rob*(Xc(n,:)-m(g,:))';
            t3=BVNcdf([t22(n,idxO(pairs.index(:,idx)),g),t11(n,idxO(pairs.index(:,idx)),g)],muc',s2c);
            %t3=mvncdf(t22(n,idxO(pairs.index(:,idx)),g),t11(n,idxO(pairs.index(:,idx)),g),muc',s2c);
            t3(t3<2.2204e-16)=2.2204e-16;
            plimix(n,g,idx)=t3*mvnpdf(Xc(n,:),m(g,:),SIGMA(idxC,idxC,g));
        end
    end
end
f=log(plimix);
%ind=(f<log(eps));
%f=(log(eps)*ind+f.*(1-ind));
f=f.*Uplusplus;
f=-sum(f(:));
end

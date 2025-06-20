function [f, plimix]=PairplikObsAll(theta,m,th_idx,idxVar,X,Q,nopairs,pairsindex)
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
n=size(Xc,1);
N=n;
Xo=X(:,idxVar==1);
O=size(Xo,2);
no_TH=ones(1,O);
for o=1:O,
no_TH(o)=sum(th_idx==o)+1;
end
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



plimix=ones(N,G,size(pairs.index,2));
%   s2c=ones(G,O);
theta.cors=corsO;
theta1=theta;
br1=threshUniv(theta1,G,O,th_idx);
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
        s2c=SIGMA(idxO(pairs.index(:,idx)),idxO(pairs.index(:,idx)),g)-SIGMA(idxC,idxO(pairs.index(:,idx)),g)'*inv(SIGMA(idxC,idxC,g))*SIGMA(idxC,idxO(pairs.index(:,idx)),g);
        [s2,r2c]=cov2corr(s2c);
        
        %for n=1:N,
            muc=mu0(g,idxO(pairs.index(:,idx)))'+SIGMA(idxC,idxO(pairs.index(:,idx)),g)'*inv(SIGMA(idxC,idxC,g))*(Xc-m(g,:))';
%             tic, t3=BVNcdf(t11(n,idxO(pairs.index(:,idx)),g),muc',s2c)-BVNcdf([t11(n,idxO(pairs.index(1,idx)),g) t22(n,idxO(pairs.index(2,idx)),g)],muc',s2c)-BVNcdf([t22(n,idxO(pairs.index(1,idx)),g) t11(n,idxO(pairs.index(2,idx)),g)],muc',s2c)+BVNcdf(t22(n,idxO(pairs.index(:,idx)),g),muc',s2c); 
%             toc
l1=(t11(:,idxO(pairs.index(:,idx)),g)-muc')./s2';
            l2=([t11(:,idxO(pairs.index(1,idx)),g) t22(:,idxO(pairs.index(2,idx)),g)]-muc')./s2';
            l3=([t22(:,idxO(pairs.index(1,idx)),g) t11(:,idxO(pairs.index(2,idx)),g)]-muc')./s2';
            l4=(t22(:,idxO(pairs.index(:,idx)),g)-muc')./s2';
            for n=1:N,
            t3=bvnlr(l1(n,1), l1(n,2), r2c(1,2) )-bvnlr(l2(n,1), l2(n,2), r2c(1,2) )-bvnlr(l3(n,1), l3(n,2), r2c(1,2) )+bvnlr(l4(n,1), l4(n,2), r2c(1,2) );
           % t3=(mvncdf(t22(n,idxO(pairs.index(:,idx)),g),t11(n,idxO(pairs.index(:,idx)),g),muc',s2c));
            t3(t3<eps)=eps;
                 t4=mvnpdf(Xc(n,:),m(g,:),SIGMA(idxC,idxC,g));
                  %t4(t4<eps)=eps;
            plimix(n,g,idx)=t3*t4;
        end
    end
end
f=log(plimix);
ind=(f<log(eps));
f=exp(log(eps)*ind+f.*(1-ind));
end
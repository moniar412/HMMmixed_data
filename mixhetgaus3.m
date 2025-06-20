function [mu,Sig,deltam,Qm,U0,Uplus,llk,bic,iter,e]=mixhetgaus3(Y,G,eps,st)
%
% ML estimte of a mix homogeneous HMM Gaussian model
% INPUT___________________________________
% Y  : data (n x T)
% G  : number of states
% st : seed of the random number generator
% OUTPUT__________________________________
% mu     : mean vector
% sig    : variances vector
% deltam : initial probabilities vector
% Gamm   : Gamm2Qm(Gamm) gives the transition probability matrices
% pm     : weights of the mixture
% Um     : posterior probabilities of the mixture
% llk    : model log-likelihood
% bic    : bic
% iter   : number of iterations
%
% 14/02/11
%
rng(st,'twister')
G2=G*G;
[n,P,T]=size(Y);
Gausy=zeros(n,T,G);
GausyExp=zeros(n,T,G);
Alpha=zeros(n,T,G);
Beta=ones(n,T,G)/G;
U=zeros(n*(T-1),G2);
Uplus=zeros(n*T,G);
onesn=ones(n,1);
onesnT=ones(n*T,1);
oP=ones(P,1);
MOmeQ=zeros(P,G*P);
MOmel=zeros(P,G);
%
Yt=Y(:,:,1);
for t=2:T,
    Yt=vertcat(Yt,Y(:,:,t));
end
U0=randU(n*T,G);
%[U0,Mmu,Sig,dif,like,bic,it]=mixhom(Yt,U0,10^-6,1);
[U0,Mmu,MOmeQ,MOmel,dif,like,awe,bic,aic,it]=mixhet(Yt,U0,0.0000001,1);
% starting point
%
deltam=sum(U0(1:n,:))'/n;
Qm=zeros(G);
for t=1:T-1
    Qm=Qm+U0(((t-1)*n+1):(t*n),:)'*U0((t*n+1):((t+1)*n),:);
end
Qm=diag(1./sum(Qm,2))*Qm;
%
dif=1;
llold=-Inf;
iter=0;
%
%
e=0;
while dif>eps
    %
    iter=iter+1;
    %if mod(iter,10)==0, disp([iter,llk]), end
    %
    % compute the Gaussian densities (n x T x G)
    for g=1:G
        Q=MOmeQ(:,(g-1)*P+(1:P));
        lam=MOmel(:,g);
        for t=1:T
            Gausy(:,t,g)=-0.5*sum(log(lam))-0.5*sum(((Y(:,:,t)-onesn*Mmu(g,:))*Q*diag(1./sqrt(lam))).^2,2);
            %ind=(Gausy(:,t,g)< -36.0437);
            %GausyExp(:,t,g)=exp( -36.0437*ind+Gausy(:,t,g).*(1-ind));
            GausyExp(:,t,g)=exp( Gausy(:,t,g));
        end
    end
    %
    % update Alpha (n x T x G) and Beta (n x T x G)
    A=squeeze(GausyExp(:,1,:)).*repmat(deltam',n,1);
    a=sum(A,2);
    um=log(a);
    Alpha(:,1,:)=repmat(1./a,1,G).*A;
    for t=2:T
        A=(squeeze(Alpha(:,t-1,:))*Qm).*squeeze(GausyExp(:,t,:));
        a=sum(A,2);
        Alpha(:,t,:)=repmat(1./a,1,G).*A;
        um=um+log(a);
        B=(Qm*squeeze(Beta(:,T-t+2,:).*GausyExp(:,T-t+2,:))')';
        b=sum(B,2);
        Beta(:,T-t+1,:)=repmat(1./b,1,G).*B;
    end
    %
    %ind=(um<log(1.0e-300));
    %um=(log(1.0e-300))*ind+um.*(1-ind);
    % update transition probabilities
    for j=1:G
        for k=1:G
            U(:,k+G*(j-1))=reshape(Alpha(:,1:T-1,k).*Beta(:,2:T,j).*GausyExp(:,2:T,j),n*(T-1),1).*Qm(k,j);
        end
    end
    U=repmat(1./sum(U,2),1,G2).*U;
    Qm=reshape(sum(U),G,G);
    Qm=repmat(1./sum(Qm,2),1,G).*Qm;
    %
    % update Uplus (n*T x G x Gm)
    Uplus=reshape(Alpha.*Beta,n*T,G);
    Uplus=repmat(1./sum(Uplus,2),1,G).*Uplus;
    %
    % update delta
    deltam=sum(Uplus(1:n,:))';
    deltam=deltam/sum(deltam);
    %
    %update mu and sig
    sw=sum(Uplus);
    mu=diag(1./sw)*Uplus'*Yt;
    %
    for g=1:G
        Sig0=Yt-onesnT*mu(g,:);
        Sig(:,:,g)=Sig0'*diag(Uplus(:,g))*Sig0/sw(g);
        [P0,L,Q]=svd(Sig(:,:,g));
        l=diag(L);
        MOmeQ(:,(g-1)*P+(1:P))=Q;
        MOmel(:,g)=min(oP*1000000,max(oP*0.000001,l));
        Sig(:,:,g)=Q*diag(MOmel(:,g))*Q';
    end
    %    
    % compute the log-likelihood
    llk=sum(um);
    dif = (llk-llold);
    llold=llk;
    %
end
%
% check number of parameters
bic = 2*llk - log(n*T)*(G-1 + G*(G-1) + G*P + G*P*(P+1)/2);
%
disp(sprintf('mixH-HMMGaus(%g): llike=%g, BIC=%g, iter=%g, dif=%g',G,llk,bic,iter,dif))
if isnan(dif)
    disp('------------------------------------------ mixhetGaus3: dif is NaN! ------------------------------------------------')
    dif
    e=1;
end
if abs(dif)>10^-4
    disp('-------------------------------------- mixhetGaus3: errore negative dif ! ------------------------------------------')
    dif
    e=1;
end

Gamm=Qm2Gamm(Qm);
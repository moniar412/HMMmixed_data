function [Y,Z,M,Sig] = genHMM_mvn0(N,T,delta,Qp,st,A,eta,Omega,L,br,th_idx,P,Q,O,idxO)
% n = # units
% T = # times
% mean = mean vector (K\times P); K = #states; P = # variables
% Sigma = covariance matrix
if isempty(N)
    N=100;
end
if isempty(T)
    T=5;
end
rng(st,'twister');
%
G = size(Qp,1);
M=zeros(G,P);
Sig=zeros(P,P,G);
for g=1:G,
    M(g,:)=(A*eta(g,:)')';
    Sig(:,:,g)=A*[Omega(:,:,g) zeros(Q,P-Q); zeros(P-Q,Q) L]*A';
end
%
Y = zeros(N,P,T);
Z = zeros(N,T);
%
% for g = 1:G
%     Sig(:,:,g) = sqrtm(Sig(:,:,g));
% end
for n = 1:N
    cp = cumsum(delta);
    for t = 1:T
        g = 1+sum(rand(1)>cp);
        Y(n,:,t) =  M(g,:) + randn(1,P)*Sig(:,:,g);
        Z(n,t) = g;
        cp = cumsum(Qp(g,:));
    end
end
b=cell(O,1);
for i=1:O,
    b{i,1}=[-Inf br(th_idx==i) Inf];
end
X1=ones(N,O,T);
for t=1:T,
    
    for i=1:O,
        
        X1(:,i,t)=ordinal(Y(:,idxO(i),t),[],[],b{i,1});
        %     ordinal(d(:,C+i),[],[],b{i,1});
    end
end
idxC=1:P;
idxC(idxO)=[];
Y(:,idxC,:)=Y(:,idxC,:);
Y(:,idxO,:)=X1;


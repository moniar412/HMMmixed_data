function ar=simula(nk)
%
%
%
ar=zeros(nk,3);
T=10; N=100;
P=8; Qt=[2 2 2 2 2]; G=3; O=4;
Q=4;
C=P-O;
idxO=[1 2 5 6];
%eta=[2 2 -1 -1 0 0 0 0;-3 -3 2 2 0  0 0 0; 5 5 5 5 0 0 0 0];
eta=[2 2 -1 -1 0 0 0 0;-1 -1 1 1 0  0 0 0; 4 4 4 4 0 0 0 0];

Omega=zeros(Q,Q,G);
Omega(:,:,1)=eye(Q);
Omega(:,:,2)=1.5*[2.5 1.2 1.2 1.2; 1.2 2.5 1.2 1.2; 1.2 1.2 2.5 1.2; 1.2 1.2 1.2 2.5];
Omega(:,:,3)=1.5*[1.5 0.75 0.75 0.75; 0.75 1.5 0.75 0.75; 0.75 0.75 1.5 0.75; 0.75 0.75 0.75 1.5];
L=eye(P-Q);
A=eye(P);
pg=[.3 .3 .4];
gamma=repmat([0 1 2 3],1,O);
th_idx=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
delta = [.3 .3 .4 ];
Qp = [.9 .05 0.05; .1 .8  .1; 0.05 .1 0.85];
no.pairs=O*(O-1)/2;
pairs.index=combnk(1:O,2)';
if O>3 pairs.index=fliplr(pairs.index); end
idxVar=[1 1  0 0 1 1 0 0];
no_TH=[5 5 5 5];
%
for r=1:nk,
    st=r+1700
    [Y,Z,M,Sig] = genHMM_mvn0(N,T,delta,Qp,st,A,eta,Omega,L,gamma,th_idx,P,Q,O,idxO);
    X=Y;
    % NO REDUCTION -   tic,[~,~,~,~,~,~,~,~,ariN]=mixhetmixedN(Y,G,O,th_idx,idxVar,idxO,no_TH,10^-4,st,Z);toc
    ar(r,2)=ariN;
    % HMM SCR
    [m,theta,Deltam,Qm,Uplus,llk,bic,iter,ari,Qm0]=mixhetmixed(Y,G,O,idxO,Q,th_idx,idxVar,no_TH,10^-3,st,Z);
%Qm
ar(r,1)=ari;
% NAIVE
[~,~,~,~,U0,~,~,~,~]=mixhetgaus3(Y,G,10^-6,st);
ar(r,3)=mrand(ftoh(U0)'*vtoh(Z(:)));
monia1=ar(r,:);
monia2=Qm(:)';
   save('Qt10n100nowell.txt','monia2','-ascii','-append')
    save('ARIt10n100nowell.txt','monia1','-ascii','-append')
end
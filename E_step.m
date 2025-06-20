function T=E_step(theta,th_idx,P,G,n)
% 
% OUTPUT__________________________________
% compute the E-step of the EM algorithm, i.e. we calculate the matrix with group membership probabilities for 
% each frquency. nrow=the total number of pi^xi_ci (it is not equal to N!!!);ncol=the total number of groups, G.
% INPUT___________________________________
% theta: a list containing the parameters of each group (by row).
% the parameter vector is composed of means,varaiances,, chol decomposition, correlations and
% the thresholds that are invariant across the groups (gamma).
% th_idx: a vector keeping track to which variable the thresholds belongs to
% n_xixj: the vector of frequencies n_xixj_cicj
% P: the number of variables 
% G:  the number of groups
% 
% to ge the variant thresholds across the groups
TH=thresh(theta,th_idx,P,G);
plik=ones(n,G);
for g=1:G,
    plik(:,g)=(PI_by_cell(TH(g,:),theta.cors(g,:),P,th_idx)).*theta.pg(g);
end
% probabilities in each row must sum to 1 (application of Bayes law)
% tic,T=bsxfun(@times, plik, 1./(sum(plik, 2)));,toc
T=repmat(1./sum(plik,2),1,G).*plik;
end


% plik1=plik.*repmat(n_xixj',1,2);
% T1=bsxfun(@times, plik1, 1./(sum(plik1, 2)))
% T_den=sum((T1.*repmat(n_xixj',1,2)),2)
% T1=T1./repmat(T_den,1,2)
% T11=T1(:,1); T11(isnan(T11))=[];
% T12=T1(:,2); T12(isnan(T12))=[];
% T1=[T11 T12]
% T2=bsxfun(@times, T1, 1./(sum(T1, 2)));
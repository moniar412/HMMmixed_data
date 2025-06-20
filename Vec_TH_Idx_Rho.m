function [th_var1_pair, th_var2_pair, rho_vector, idx_th_var1_pair, idx_th_var2_pair, idx_var1_pair, idx_var2_pair, idx_extended_pairs, last_th_var1_pair, last_th_var2_pair, th1_dnorm_in_der_pi_to_gamma_xi, th2_dnorm_in_der_pi_to_gamma_xj, th2_minus_last_th1, th1_minus_last_th2]=Vec_TH_Idx_Rho(P, TH, th_idx, cors)
%
% OUTPUT__________________________________
% it is a list of various vectors, either vectors keeping track of
% indices or vectors with thresholds and rho's duplicated properly, all needed
% for the computation of expected probabilities and derivatives.
% all duplications are done as follows: within each pair of variables,xixj, if
% for example we want to duplicate their thresholds, tau^xi_ci and tau^xj_cj, 
% then index ci runs faster than cj (for each cj we take all different tau^ci's,
% and then we proceed to next cj and do the same).
% all pairs xixj are taken with index j running faster than i
% each variable may have a different number of categories 
% INPUT___________________________________
% no.x is the number of ordinal variables
% TH is a vector containing the thresholds of all variables (th_var1,th_var1,...,
% th_varP,th_varP). TH does not contain tau_0=-Inf and tau_last=Inf for all
% variables.
% th.idx is a vector keeping track to which variable the thresholds in TH 
% belongs to
% cors is the vector of all correlations (rho_xixj) where j runs faster than i  
% 
% 
% number of thresholds for each variable

for i=1:P,
    no.TH(i)=length(TH(th_idx==i));
end
pairs.index=combnk(1:P,2)';
if P>3 pairs.index=fliplr(pairs.index); end
no.pairs=size(pairs.index,2);
% 
% 
% all the possible combinations of the thresholds of any two variables
% the combinations(-Inf,-Inf),(-Inf,Inf),(Inf,-Inf) and (Inf,Inf) are
% NOT included
%  
 th_var1_pair=cell(no.pairs,1);
 th_var2_pair=cell(no.pairs,1);
 rho_vector=cell(no.pairs,1);
%
%  vectors containing the indices of all thresholds (from -Inf to Inf)
%  the voectors are longer than the vectors TH.var1.pair & TH.var2.pair
% 
 idx_th_var1_pair =cell(no.pairs,1);
 idx_th_var2_pair =cell(no.pairs,1);
%  
%  vectors keeping track to which variable the thresholds in TH belongs to
%  the two vectors are of same length as that of idx.TH.var1.pair & idx.TH.var2.pair
%  
 idx_var1_pair =cell(no.pairs,1);
 idx_var2_pair =cell(no.pairs,1);
%  
%  
%  the index of the pair for each pair of variables - example: pair of
%  variables 1-2 has index 1, variables 1-3 has index 2...
%  
idx_extended_pairs =cell(no.pairs,1);
% 
%  th1.for.dnorm.in.der.pi.to.tau.xi and
%  th2.for.dnorm.in.der.pi.to.tau.xj give the thresholds of almost
%  all variables appropriately duplicated so that the vectors can be used
%  as input in dnorm() to compute der.pi.xixj.to.tau.xi and
%  der.pi.xixj.to.tau.xj.
%  th1.for.dnorm.in.der.pi.to.tau.xi does not contain the thresholds of
%  the last variable and th2.for.dnorm.in.der.pi.to.tau.xj those of
%  the first variable
% 
th1_dnorm_in_der_pi_to_gamma_xi=cell(no.pairs,1); 
th2_dnorm_in_der_pi_to_gamma_xj =cell(no.pairs,1);

for i=1:no.pairs,
   one.th.var1.pair=TH(th_idx==pairs.index(1,i));
   one.th.var2.pair=TH(th_idx==pairs.index(2,i));

%    recall that the first (-Inf) and last (Inf) thresholds are not included
%    thus no.th.var1.pair is equal to the number of categories of var1 minus 1
%    similarly for no.th.var2.pair

   no.th.var1.pair=no.TH(pairs.index(1,i));
   no.th.var2.pair=no.TH(pairs.index(2,i));
   th_var1_pair{i,1}=repmat(one.th.var1.pair,1,no.th.var2.pair );
   th_var2_pair{i,1}=reshape((repmat(one.th.var2.pair',1,no.th.var1.pair))',1,no.th.var1.pair*length(one.th.var2.pair));   
   rho_vector{i,1}=repmat(cors(i),1,length(th_var1_pair{i,1}));
   idx_th_var1_pair{i,1}=repmat(0:(no.th.var1.pair+1),1,(no.th.var2.pair+2));
   idx_th_var2_pair{i,1}= reshape((repmat((0:(no.th.var2.pair+1))',1,(no.th.var1.pair+2)))',1,length(0:(no.th.var2.pair+1))*(no.th.var1.pair+2));
   length_vec=length(idx_th_var1_pair{i,1});
   idx_var1_pair{i,1}=repmat(pairs.index(1,i),1,length_vec);
   idx_var2_pair{i,1}=repmat(pairs.index(2,i),1,length_vec);
   idx_extended_pairs{i,1}=repmat(i,1,length_vec);
   th1_dnorm_in_der_pi_to_gamma_xi{i,1}=repmat(one.th.var1.pair,1,(no.th.var2.pair+1));
   th2_dnorm_in_der_pi_to_gamma_xj{i,1}=reshape((repmat(one.th.var2.pair',1,(no.th.var1.pair+1)))',1,length(one.th.var2.pair)*(no.th.var1.pair+1));
end
%  fill the vectors
%is there a different way to change the dimensions?
th_var1_pair=reshape((cell2mat(th_var1_pair'))',1,size(cell2mat(th_var1_pair'),1)*size(cell2mat(th_var1_pair'),2)); 
th_var2_pair=reshape((cell2mat(th_var2_pair'))',1,size(cell2mat(th_var2_pair'),1)*size(cell2mat(th_var2_pair'),2)); 
rho_vector=reshape((cell2mat(rho_vector'))',1,size(cell2mat(rho_vector'),1)*size(cell2mat(rho_vector'),2));
idx_th_var1_pair=reshape((cell2mat(idx_th_var1_pair'))',1,size(cell2mat(idx_th_var1_pair'),1)*size(cell2mat(idx_th_var1_pair'),2));
idx_th_var2_pair=reshape((cell2mat(idx_th_var2_pair'))',1,size(cell2mat(idx_th_var2_pair'),1)*size(cell2mat(idx_th_var2_pair'),2));
idx_var1_pair=reshape((cell2mat(idx_var1_pair'))',1,size(cell2mat(idx_var1_pair'),1)*size(cell2mat(idx_var1_pair'),2));
idx_var2_pair=reshape((cell2mat(idx_var2_pair'))',1,size(cell2mat(idx_var2_pair'),1)*size(cell2mat(idx_var2_pair'),2));
idx_extended_pairs=reshape((cell2mat(idx_extended_pairs'))',1,size(cell2mat(idx_extended_pairs'),1)*size(cell2mat(idx_extended_pairs'),2));
th1_dnorm_in_der_pi_to_gamma_xi=reshape((cell2mat(th1_dnorm_in_der_pi_to_gamma_xi'))',1,size(cell2mat(th1_dnorm_in_der_pi_to_gamma_xi'),1)*size(cell2mat(th1_dnorm_in_der_pi_to_gamma_xi'),2));
th2_dnorm_in_der_pi_to_gamma_xj=reshape((cell2mat(th2_dnorm_in_der_pi_to_gamma_xj'))',1,size(cell2mat(th2_dnorm_in_der_pi_to_gamma_xj'),1)*size(cell2mat(th2_dnorm_in_der_pi_to_gamma_xj'),2));
% 
% vectors needed to compute the expected probabilities.
% when the first variable of the pair has gamma_last=Inf and the second
% a finite threshold then compute "pnorm(th2.minus.last.th1)".
% when the second variable of the pair has gamma_last=Inf and the first
% a finite threshold then compute "pnorm(th1.minus.last.th2)".
% 
 th2_minus_last_th1=cell((P-1),1);
 th1_minus_last_th2=cell((P-1),1);
 for i=1:(P-1),
     th2_minus_last_th1{i}= TH(ismember(th_idx,((i+1):P)));
     th1_minus_last_th2{i}=repmat(TH(th_idx==i),1,(P-i));
  end
%
% fill the vectors
 th2_minus_last_th1 =cell2mat(th2_minus_last_th1');
 th1_minus_last_th2 =cell2mat(th1_minus_last_th2');
%
% the vector includes T/F showing which elements of idx.TH.var1.pair
% correspond to the last thresholds of variables.
% the length is the same as that of idx.TH.var1.pair 

last_th_var1_pair=idx_var1_pair==1 & idx_th_var1_pair==(no.TH(1)+1);
%
% we consider (no.x-1) variables because the first one is not included 
% (j runs faster than i)

for i=2:(P-1),
   condition=idx_var1_pair==i & idx_th_var1_pair==(no.TH(i)+1);
   last_th_var1_pair= (last_th_var1_pair | condition);
end

%  the vector includes T/F showing which elements of idx.TH.var2.pair
%  correspond to the last thresholds of variables.
%  since in pairs xi.xj, where j runs faster than i, the first variable is not
%  included in the comlumn of xj's, then we start with the second variable.
%  the length is the same as that of idx.TH.var2.pair 
last_th_var2_pair=idx_var2_pair==2 & idx_th_var2_pair==(no.TH(2)+1);
 for i=3:P,
     condition=idx_var2_pair==i & idx_th_var2_pair==(no.TH(i)+1);
    last_th_var2_pair = (last_th_var2_pair | condition);
 end

end







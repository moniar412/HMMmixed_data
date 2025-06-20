function PWexp = PWexpected(th_var1_pair, th_var2_pair, rho_vector, idx_th_var1_pair, idx_th_var2_pair, idx_var1_pair, idx_var2_pair, idx_extended_pairs, last_th_var1_pair, last_th_var2_pair, th1_dnorm_in_der_pi_to_gamma_xi, th2_dnorm_in_der_pi_to_gamma_xj, th2_minus_last_th1, th1_minus_last_th2)
%
% OUTPUT__________________________________
%pi^xixj_cicj is the expected probability for the pair of variables xixj and categories cicj.
%The index ci runs the fastest of all, followed by cj, then by j and lastly by i.
% INPUT___________________________________
% the output of  Vec.TH.Idx.Rho
% 
prob.vec=repmat(nan,1,length(idx_th_var1_pair));
prob.vec(idx_th_var1_pair==0 | idx_th_var2_pair==0)=0; 
prob.vec(last_th_var1_pair & last_th_var2_pair)=1;
prob.vec(last_th_var1_pair & idx_th_var2_pair~=0 & ~last_th_var2_pair)=cdf('Normal',th2_minus_last_th1,0,1);
prob.vec(last_th_var2_pair & idx_th_var1_pair~=0 & ~last_th_var1_pair)=cdf('Normal',th1_minus_last_th2,0,1);
%Can I avoid the for loop???
a=zeros(length(th_var1_pair),1);
for i=1:length(th_var1_pair)
    a(i)=bvnlr(th_var1_pair(i),th_var2_pair(i),rho_vector(i));
end
prob.vec(isnan(prob.vec))=a;
cum1=prob.vec(idx_th_var1_pair~=0 & idx_th_var2_pair ~=0);
cum2=prob.vec(idx_th_var1_pair~=0 & ~last_th_var2_pair);
cum3=prob.vec(idx_th_var2_pair ~=0 & ~last_th_var1_pair);
cum4=prob.vec(~last_th_var1_pair & ~ last_th_var2_pair);
PWexp=cum1-cum2-cum3+cum4 ;
%
% to avoid Nan/-Inf: log(.Machine$double.eps)=-36.04365
% all elements should be strictly positive
%
PWexp(PWexp<eps)=eps;
end



 



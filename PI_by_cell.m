function PI = PI_by_cell(TH,cors,P,th_idx)
%
% OUTPUT__________________________________
% the expected probability for the pair of variables xixj and thecategories cicj.
% index ci runs the fastest of all,followed by cj,then by j and lastly by i.
% INPUT___________________________________
% TH: the vector of thresholds Tau (varying across the groups)
% cors: the correlations vector:rho_xixj where j runs faster than i.
% th_idx: the vector keeping track to which variable the thresholds belongs to
% P: the number of ordinal variables 
% 
[th_var1_pair, th_var2_pair, rho_vector, idx_th_var1_pair, idx_th_var2_pair, idx_var1_pair, idx_var2_pair, idx_extended_pairs, last_th_var1_pair, last_th_var2_pair, th1_dnorm_in_der_pi_to_gamma_xi, th2_dnorm_in_der_pi_to_gamma_xj, th2_minus_last_th1, th1_minus_last_th2]=Vec_TH_Idx_Rho(P, TH, th_idx, cors);
PI=PWexpected(th_var1_pair, th_var2_pair, rho_vector, idx_th_var1_pair, idx_th_var2_pair, idx_var1_pair, idx_var2_pair, idx_extended_pairs, last_th_var1_pair, last_th_var2_pair, th1_dnorm_in_der_pi_to_gamma_xi, th2_dnorm_in_der_pi_to_gamma_xj, th2_minus_last_th1, th1_minus_last_th2);

end

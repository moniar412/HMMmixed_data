function TH=thresh(theta,th_idx,P,G)
% 
% OUTPUT__________________________________
% a matrix containing the thresholds (Tau) varying across the groups.
% the rows are the thresholds in each group. The order in each row is th_var1,th_var1,...,th_varP,th_varP.
% INPUT___________________________________
% theta: a list of parameters(gamma,means, varaiances, chol decompositiona and correlations) in each group
% th_idx: a vector keeping track to which variable the thresholds belongs to
% G: the number of gropus
% P: the number of ordinal variables
% 
% to compute the unconditional mean and covariance matrix
% mug=theta.pg*theta.mu;
% var.W=theta.pg*theta.sigma;
% % to compute var.B
%  m=reshape((theta.mu'-repmat(mug,2,1)'),P,G);
% % m=theta.mu'-repmat(mug,2,1)';
% %matrix of variances: PxG
% var.B=sum((m.^2).*repmat(theta.pg,P,1),2);
% % to compute var.W
% var.G=var.W+var.B';
% %IF THE VARIABLES HAVE THE SAME NUMBER OF CATEGORIES...how can we
% %generalize it???
% br=mug(th_idx)+theta.gamma.*sqrt(var.G(th_idx));
br=theta.gamma;
br1=repmat(br,G,1);

TH=(br1- theta.mu(:,th_idx))./sqrt(theta.sigma(:,th_idx));
end
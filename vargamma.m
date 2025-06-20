function [gamma1_P,gamma1_Pext1,gamma1_Pext2]=vargamma(gamma,th_idx,P)
%
% it is a general function needed when the variables have both different
% number of categories and different threshold values.
% the output
%
gamma1_P=cell(P,1);
gamma1_Pext1=cell(P,1);
gamma1_Pext2=cell(P,1);
for p=1:P,
    gamma1_P{p}=sort(gamma(th_idx==p));
    gamma1_Pext1{p}=[gamma1_P{p} +Inf];
    gamma1_Pext2{p}=[-Inf gamma1_P{p}];
end
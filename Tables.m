function tables=Tables(X,P,no_TH)
%
% OUTPUT__________________________________
% tables: two-way frequency table
% INPUT___________________________________
% X: the data matrix
% P: the number of ordinal variables
%
% pairs.index=combntns(1:P,2)';
pairs.index=combnk(1:P,2)';
if P>3 pairs.index=fliplr(pairs.index); end
no.pairs=size(pairs.index,2);
v=pairs.index(1,1:no.pairs);
u=pairs.index(2,1:no.pairs);
tables=cell(no.pairs,1);
for i=1:no.pairs,
    tables{i,1}=crosstab(X(:,v(i)),X(:,u(i)));
    if (size(tables{i},1)~=no_TH(v(i))),  tables{i}=[tables{i};zeros(no_TH(v(i))-size(tables{i},1),size(tables{i},2))]; 
        end
    if (size(tables{i},2)~=no_TH(u(i))),  tables{i}=[tables{i},zeros(size(tables{i},1),no_TH(u(i))-size(tables{i},2))]; end   
end
end

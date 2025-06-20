function n_xixj=frequency1(tables,P)
no.pairs=P*(P-1)/2;
n_xixj=cell(no.pairs,1);
for i=1:no.pairs,
n_xixj{i}=reshape(tables{i},1,[]);
end
% n_xixj=reshape(cell2mat(n_xixj)',1,[]);
n_xixj=reshape(cell2mat(n_xixj'),1,[]);
end
function [e,n]=expectationvalue(mps,hset)

[M,N]=size(hset); 
d=size(mps{1},3);

% expectation value 
e=0;
for m=1:M
    em=1;
    for j=N:-1:1
        h=hset{m,j};
        h=reshape(h,[1,1,d,d]); 
        em=updateCright(em,mps{j},h,mps{j});
    end
    e=e+em; 
end

% norm
n=1;
X=eye(d); X=reshape(X,[1,1,d,d]); 
for j=N:-1:1
    n=updateCright(n,mps{j},X,mps{j}); 
end

e=e/n;
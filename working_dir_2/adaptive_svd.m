function [U,S,V,listn,j] = adaptive_svd(mps,mpo,site,p,iter,r,epsilon)
%
%random SVD with adaptive range finder
%
[D,~,~] = size(mps{site});
[x,~,~,~] = size(mpo{site});
%
%find the basis Q which satisfies a given error
%
[q,listn,j] = adaptive_rand(mps,mpo,site,r,epsilon,iter);

l= size(q,2);
B = zeros(l,D*x*D);
for row = 1:l
    aux = mpoxrand(mpo,mps,site,'left',q(:,row)','normal');
    B(row,:) = reshape(aux,[D*x*D,1]);
end
[U,S,V] = svd2(B);
U = U(:,1:p);
S = S(1:p,1:p);
V = V(1:p,:);
V = S*V;
U = q*U;
U = reshape(U,[D,x,D,p]);
V = reshape(V,[p,D,x,D]);


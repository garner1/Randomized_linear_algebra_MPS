function [U,V] = random_svd_h(mps,mpo,site,p)
%
%Perform random svd for the effective Hamiltonian
%
[D,~,~] = size(mps{site});
Y = zeros(D^2*5,p);
for col = 1:p
    p0 = randn(D,5,D)+1i*randn(D,5,D);
    aux = mpoxrand(mpo,mps,site,'right',p0,'normal');
    Y(:,col) = reshape(aux,[D*5*D,1]);
end
[q,~] = qr(Y,0);

B = zeros(p,D^2*5);
for row = 1:p
    aux = mpoxrand(mpo,mps,site,'left',q(:,row)','normal');
    B(row,:) = reshape(aux,[D*5*D,1]);
end
[U,Sh,V] = svd2(B);
V = Sh*V; 
U = q*U;
U = reshape(U,[D,5,D,p]);
V = reshape(V,[p,D,5,D]);

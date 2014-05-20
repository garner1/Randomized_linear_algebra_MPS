function [U,S] = random_eig(mps,mpo,site,p,iter,l)
%
%Perform random svd for the normalization
%
[D,~,~] = size(mps{site});
[x,~,~,~] = size(mpo{site});
Y = zeros(D*x*D,p);
%
%find the basis Q which satisfies a given error
%
% r = 10;
% [q,~] = adaptive_rand(mps,mpo,site,r,1e-10,iter);

for col = 1:l
    p0 = randn(D,x,D)+1i*randn(D,x,D);
    aux = mpoxrand(mpo,mps,site,'right',p0,'normal');
    Y(:,col) = reshape(aux,[D*x*D,1]);
end
[q,~] = qr(Y,0);
for ind = 1:iter
    for col = 1:l
        p0 = reshape(q(:,col),[D,x,D]);
        aux = mpoxrand(mpo,mps,site,'right',p0,'dagger');
        Y(:,col) = reshape(aux,[D*x*D,1]);
    end
    [q,~] = qr(Y,0);
    for col = 1:l
        p0 = reshape(q(:,col),[D,x,D]);
        aux = mpoxrand(mpo,mps,site,'right',p0,'normal');
        Y(:,col) = reshape(aux,[D*x*D,1]);
    end
    [q,~] = qr(Y,0);
end
% l= size(q,2);
B = zeros(l);
for row = 1:l
    for col = 1:l
        aux = mpoxrand(mpo,mps,site,'left',q(:,row)','normal');
        aux = reshape(aux,[1,D*x*D]);
        B(row,col) = aux*q(:,col);
    end
end
[U,S] = eig(B);
U = U(:,1:p);
S = S(1:p,1:p);
U = q*U;
U = reshape(U,[D,x,D,p]);


function [U,S,V] = random_svd(mps,mpo,site,p,iter,l)
%
%Perform random svd for the normalization
% site=site to be excluded from the network
% p=number of needed singular values
% iter=number of iterations to increase accuracy
% l=number of random states generated
%
[D,~,~] = size(mps{site});
[x,~,~,~] = size(mpo{site});
Y = zeros(D*x*D,p);

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
B = zeros(l,D*x*D);
for row = 1:size(q,2)
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


function [U,V] = randomEnv2(mps,mpo,p,iter,l)
%
% p=number of needed singular values
% iter=number of iterations to increase accuracy
% l=number of random states generated
%
[D,~,~] = size(mps{1});
[x,~,~,~] = size(mpo{1});
Y = zeros(D*x*D,p);
%
%compute Q
%
for col = 1:2*l
    Y(:,col) = reshape(mpoxrand(mpo,mps,'right',randn(D,x,D)+1i*randn(D,x,D),'normal'),[D*x*D,1]);
end
[q,~] = qr(Y,0);
%
%compute B
%
B = zeros(2*l,D*x*D);
for row = 1:size(q,2)
    B(row,:) = reshape(mpoxrand(mpo,mps,'left',q(:,row)','normal'),[D*x*D,1]);
end
%
%svd of B
%
[Omega,~,~] = svd2(B);
Omega = Omega(:,1:l);
Omega = q*Omega;
%
%second part
%
l = size(Omega,2);
for col = 1:l
    Y(:,col) = reshape(mpoxrand(mpo,mps,'right',reshape(Omega(:,col),[D,x,D]),'normal'),[D*x*D,1]);
end
[q,~] = qr(Y,0);
for ind = 1:iter
    for col = 1:l
        Y(:,col) = reshape(mpoxrand(mpo,mps,'right',reshape(q(:,col),[D,x,D]),'dagger'),[D*x*D,1]);
    end
    [q,~] = qr(Y,0);
    for col = 1:l
        Y(:,col) = reshape(mpoxrand(mpo,mps,'right',reshape(q(:,col),[D,x,D]),'normal'),[D*x*D,1]);
    end
    [q,~] = qr(Y,0);
end
B = zeros(l,D*x*D);
for row = 1:size(q,2)
    B(row,:) = reshape(mpoxrand(mpo,mps,'left',q(:,row)','normal'),[D*x*D,1]);
end
yy = q(:,1:p);
bb = B(1:p,:);
U = reshape(yy,[D,x,D,p]);
V = reshape(bb,[p,D,x,D]);
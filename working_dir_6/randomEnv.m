function [U,V] = randomEnv(mps,mpo,p,iter,l)
%
% p=number of needed singular values
% iter=number of iterations to increase accuracy
% l=number of random states generated
%
[D,~,~] = size(mps{1});
[x,~,~,~] = size(mpo{1});
Y = zeros(D*x*D,p);

for col = 1:l
    vec0 = randn(D,x,D)+1i*randn(D,x,D);
    aux = mpoxrand(mpo,mps,'right',vec0,'normal');
    Y(:,col) = reshape(aux,[D*x*D,1]);
end
[q,~] = qr(Y,0);
for ind = 1:iter
    for col = 1:l
        vec0 = reshape(q(:,col),[D,x,D]);
        aux = mpoxrand(mpo,mps,'right',vec0,'dagger');
        Y(:,col) = reshape(aux,[D*x*D,1]);
    end
    [q,~] = qr(Y,0);
    for col = 1:l
        vec0 = reshape(q(:,col),[D,x,D]);
        aux = mpoxrand(mpo,mps,'right',vec0,'normal');
        Y(:,col) = reshape(aux,[D*x*D,1]);
    end
    [q,~] = qr(Y,0);
end
B = zeros(l,D*x*D);
for row = 1:size(q,2)
    aux = mpoxrand(mpo,mps,'left',q(:,row)','normal');
    B(row,:) = reshape(aux,[D*x*D,1]);
end
yy = q(:,1:p);
bb = B(1:p,:);
U = reshape(yy,[D,x,D,p]);
V = reshape(bb,[p,D,x,D]);
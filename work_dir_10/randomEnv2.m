function [U,V] = randomEnv2(mps,mpo,p,iter,l)
%
% p=number of needed singular values
% iter=number of iterations to increase accuracy
% l=number of random states generated
%
[D,~,~] = size(mps{1});
[x,~,~,~] = size(mpo{1});
%
%compute Q
%
Y = reshape(cell2mat(cellfun(@(x) mpoxrand(mpo,mps,'right',x,'normal'),num2cell(randn(D,x,D,2*l)+1i*randn(D,x,D,2*l),[1 2 3]),'UniformOutput',false)),[D*x*D,2*l]);
[q,~] = qr(Y,0);
%
%compute B
%
aux1 = cellfun(@(x) mpoxrand(mpo,mps,'left',x,'normal'),num2cell(conj(q),1)','UniformOutput',false);
aux2 = cellfun(@(xx) reshape(xx,[1,D*x*D]),aux1,'UniformOutput',false);
B = reshape(cell2mat(aux2),[size(q,2),D*x*D]);
%
%svd of B
%
[UU,~,~] = svd2(B);
Omega = q*UU(:,1:l);
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
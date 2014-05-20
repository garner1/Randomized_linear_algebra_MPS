function [U,V] = randomEnv(mps,mpo,p,iter,l)
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
Y = reshape(cell2mat(cellfun(@(x) mpoxrand(mpo,mps,'right',x,'normal'),num2cell(randn(D,x,D,l)+1i*randn(D,x,D,l),[1 2 3]),'UniformOutput',false)),[D*x*D,l]);
[q,~] = qr(Y,0);
for ind = 1:iter
    Y = reshape(cell2mat(cellfun(@(x) mpoxrand(mpo,mps,'right',x,'dagger'),num2cell(reshape(q,[D,x,D,l]),[1 2 3]),'UniformOutput',false)),[D*x*D,l]);
    [q,~] = qr(Y,0);
    Y = reshape(cell2mat(cellfun(@(x) mpoxrand(mpo,mps,'right',x,'normal'),num2cell(reshape(q,[D,x,D,l]),[1 2 3]),'UniformOutput',false)),[D*x*D,l]);
    [q,~] = qr(Y,0);
end
aux1 = cellfun(@(x) mpoxrand(mpo,mps,'left',x,'normal'),num2cell(conj(q),1)','UniformOutput',false);
aux2 = cellfun(@(xx) reshape(xx,[1,D*x*D]),aux1,'UniformOutput',false);
B = reshape(cell2mat(aux2),[size(q,2),D*x*D]);

yy = q(:,1:p);
bb = B(1:p,:);
U = reshape(yy,[D,x,D,p]);
V = reshape(bb,[p,D,x,D]);
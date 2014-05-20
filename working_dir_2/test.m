% clear all;
% N = 2^7;
% D = 2^4;
% mps = cell(1,N); 
% % prepare the fhe other A matrices except for the last one
% for ind = 1:N
%     U = randU(2*D);
%     U = U(:,1:D);
%     U = reshape(U,[2,D,D]);
%     mps{ind} = permute(U,[2,3,1]);
% end

mpox = mpo_id(N,2); 
site=1;
N = size(mpox,2);
[x,~,~,~] = size(mpox{1});
[m,~,~] = size(mps{1});
T = contracttensors(mpox{N},4,4,mps{N},3,3);
T = contracttensors(T,5,[3],conj(mps{N}),3,[3]);
T0 = permute(T,[5,1,3,6,2,4]);
for ind = N-1:-1:2    
    T = contracttensors(mpox{ind},4,4,mps{ind},3,3);
    T = contracttensors(T,5,[3],conj(mps{ind}),3,[3]);
    T = permute(T,[5,1,3,6,2,4]);
    T0 = contracttensors(T,6,[4,5,6],T0,6,[1,2,3]);
end
p = reshape(T0,[m*x*m,m*x*m]);
pp = reshape(permute(reshape(p,[m,m,m,m]),[1,3,2,4]),[m*m,m*m]);

[U,S,V] = svd(p);
disp(diag(S(1:D,1:D)))
disp(S(2,2))

x=contracttensors(conj(mps{site}),3,3,mps{site},3,3);
x=permute(x,[1,3,2,4]);
x=reshape(x,[D^2,D^2]);
disp(trace(p*x)); 

% j = randi(site);
% x=contracttensors(conj(mps{j}),3,[1,3],mps{j},3,[1,3]);
% x=reshape(x,[D,D]);
% disp(norm(x-eye(size(x))))

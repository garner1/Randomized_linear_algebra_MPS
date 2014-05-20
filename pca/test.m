n = 2^10;

A = rand_vec(n); A = A*A';
% A = rand(n);

k = 5;
l = k+5;
its = 5;
tic
[U,S,V] = pca_my(A,k,its,l);
S = diag(S);
toc

tic
[UU,SS,VV] = svd(A);
UU = UU(:,1:k);
SS = diag(SS);
SS = SS(1:k);
toc

disp([S(1:k),SS(1:k)])

disp([U(1:4,1),UU(1:4,1)])

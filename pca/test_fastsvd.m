i = 7;
M= 2^(i)*5;
N= 2^(i+1)*5;
d = 2;
its = 3;

A = randn(M,N,d); 
A = reshape(permute(A,[1,3,2]),[d*M,N]);

tic
[U0,S0,V0]=svd2(A);
toc
diffsnorm(A,U0(:,1:2^i),S0(1:2^i,1:2^i),V0(1:2^i,:)')

% %
% %   THIS IS THE WORSE
% %
% tic
% [U1,S1,V1]=svds(A,2^i);
% toc
% diffsnorm(A,U1,S1,V1)

tic
[U2,S2,V2]=pca2(A,2^i,its);
toc
diffsnorm(A,U2,S2,V2)

tic
[U3,S3,V3]=pca_my(A,2^i,its);
toc
diffsnorm(A,U3,S3,V3)


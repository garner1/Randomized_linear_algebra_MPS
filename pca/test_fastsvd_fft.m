i = 8;
D1= 2^(i+1)*5;
D2= 2^(i)*5;
d = 2;
its = 0;

A = randn(D1,D2,d);
A = reshape(permute(A,[1,3,2]),[d*D1,D2]);

k = min(size(A));
l = k+0;

% [U1,S1,V1]=pca_my(A,k,its,l);


[Q2,R2]=svd_fft(A,l);

% disp([diffsnorm(A,U1,S1,V1),diffsnorm(A,U3,S3,V3)])



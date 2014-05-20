clear all;

n = 4000;
k = 1;
its = 20;
l = 5;

A = randn(n);

H = A+A';

tic
[V1,D1] = eigs(H,k,'sa');
toc

tic
[V2,D2] = fast_eig(H,k,its,l);
toc

%
% DETERMINE THE POSITION OF THE APPROX MIN EIGENVALUE
%
pos = find(diag(D2)==min(diag(D2)));

% disp(pos)
% disp([V1(1:10,1),V2(1:10,pos)])
disp([D1,min(diag(D2))])
% disp(norm(V1-V2(:,pos)))

%
% ONE CAN HAVE A GOOD APPROX USING LARGE its AND SMALL l, FASTER THAN
% LANCZOS
%



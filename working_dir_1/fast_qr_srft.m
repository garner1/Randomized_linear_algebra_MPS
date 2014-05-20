function  [Q,R] = fast_qr_srft(A,k,l)

%
% Retrieve the dimensions of A.
%
[~,n] = size(A);

%
% Apply A to a random matrix, obtaining H.
%
colindices = randperm(n);
colindices = colindices(1:l);

d = 1-2*(rand(1,n) > .5);

%
% fast A*D
%
Y = fft(bsxfun(@times,A,d)'); 
F = Y(colindices, :)';
%
% Form a matrix Q whose columns constitute an orthonormal basis
% for the columns of F.
%
[Q,~] = qr(F,0);
Q = Q(:,1:k);
R = Q'*A;
clear F;


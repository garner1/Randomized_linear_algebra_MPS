clear;
m = 10^4;
n = 10^3;
k = 10;
l = 2*k;


A = randn(m,n);

Yaux = ADmult(A);

Yaux = applyfft(Yaux);

clist = randi(n,1,l);
Y = zeros(m,l);
for ind = 1:l
    Y(:,ind) = Yaux(:,clist(ind));
end

[Q,R] = qr(Y,0);
%
% SVD Q'*A to obtain approximations to the singular values
% and right singular vectors of A; adjust the left singular vectors
% of Q'*A to approximate the left singular vectors of A.
%
[U2,S,V] = svd(Q'*A,'econ');
U = Q*U2;

clear U2;

%
% Retain only the leftmost k columns of U, the leftmost k columns of V,
% and the uppermost leftmost k x k block of S.
%
U = U(:,1:k);
V = V(:,1:k);
S = S(1:k,1:k);

disp([norm(A-Q*Q'*A),norm(A-U*S*V')])



%    U = randU(n) generates a random  unitary matrix, of dimension n,
%    distributed uniformly according to the Haar measure

function U = randU(n)

X = (randn(n) + 1i*randn(n))/sqrt(2);
[Q,R] = qr(X);
R = diag(diag(R)./abs(diag(R)));
U = Q*R;
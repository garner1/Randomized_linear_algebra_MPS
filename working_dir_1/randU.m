function U = randU(n)

%    U = RANDU(n) generates a random  unitary matrix with n rows,
%    distributed uniformly according to the Haar measure.

X = (randn(n) + 1i*randn(n))/sqrt(2);
[Q,R] = qr(X);
R = diag(diag(R)./abs(diag(R)));
U = Q*R;
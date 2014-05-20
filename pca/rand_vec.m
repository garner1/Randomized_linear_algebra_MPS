%    U = randU(n) generates a random  unitary matrix, of dimension n,
%    distributed uniformly according to the Haar measure

function X = rand_vec(n)

X = (randn(n,1) + 1i*randn(n,1));
X = X/norm(X);
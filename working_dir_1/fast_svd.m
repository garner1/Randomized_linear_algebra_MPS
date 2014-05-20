function  [U,S,V] = fast_svd(A,k,its,l)

% [Q,R] = svd_fast(A,k,its,l)
% [U,S,V] = svd_fast(A,k,its,l)

%
% Retrieve the dimensions of A.
%
[m,n] = size(A);
%
% Apply A to a random matrix, obtaining H.
%
if(isreal(A))
H = A*randn(n,l);
end
if(~isreal(A))
H = A*( randn(n,l) + 1i*randn(n,l) );
end
%
% Initialize F to its final size and fill its leftmost block with H.
%
F = zeros(m,(its+1)*l);
F(1:m, 1:l) = H;
%
% Apply A*A' to H a total of its times,
% augmenting F with the new H each time.
%
for it = 1:its
H = (H'*A)';
H = A*H;
F(1:m, (1+it*l):((it+1)*l)) = H;
end

clear H;
%
% Form a matrix Q whose columns constitute an orthonormal basis
% for the columns of F.
%
[Q,~] = qr(F,0);  
%   Q = Q(:,1:k);
%   R = R(1:k,:);

clear F;

%
% SVD Q'*A to obtain approximations to the singular values
% and right singular vectors of A; adjust the left singular vectors
% of Q'*A to approximate the left singular vectors of A.
%
[U2,S,V] = svd(Q'*A,'econ');
U = Q*U2;

clear Q U2;
%
% Retain only the leftmost k columns of U, the leftmost k columns of V,
% and the uppermost leftmost k x k block of S.
%
U = U(:,1:k);
V = V(:,1:k); V = V';
S = S(1:k,1:k);

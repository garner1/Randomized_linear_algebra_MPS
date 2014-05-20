function  [Q,R] = qr_fast(A,k,its,l)
%
%   Reference:
%   Nathan Halko, Per-Gunnar Martinsson, and Joel Tropp,
%   Finding structure with randomness: Stochastic algorithms
%   for constructing approximate matrix decompositions,
%   arXiv:0909.4061 [math.NA; math.PR], 2009
%   (available at http://arxiv.org).
%
%

%
% Set the inputs k, its, and l to default values, if necessary.
%
if(nargin == 1)
  k = 6;
  its = 2;
  l = k+2;
end
if(nargin == 2)
  its = 2;
  l = k+2;
end
if(nargin == 3)
  l = k+2;
end
%
% Check the first input argument.
%
if(~isfloat(A))
  error('MATLAB:pca:In1NotFloat',...
        'Input 1 must be a floating-point matrix.')
end
if(isempty(A))
  error('MATLAB:pca:In1Empty',...
        'Input 1 must not be empty.')
end
%
% Retrieve the dimensions of A.
%
[m,n] = size(A);
% %
% % SVD A directly if (its+1)*l >= m/1.25 or (its+1)*l >= n/1.25.
% %
% if(((its+1)*l >= m/1.25) || ((its+1)*l >= n/1.25))
% % disp('slow svd')
% 
%   if(~issparse(A))
%     [U,S,V] = svd2(A);
%   end
% 
%   if(issparse(A))
%     [U,S,V] = svd2(full(A));
%   end
% %
% % Retain only the leftmost k columns of U, the leftmost k columns of V,
% % and the uppermost leftmost k x k block of S.
% %
%   U = U(:,1:k);
%   V = V(1:k,:);
%   S = S(1:k,1:k);
% 
%   return
% 
% end
% 

% if(m >= n)
% disp('fast svd')
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
  Q = Q(:,1:k);
  R = Q'*A;
  clear F;
% end

% %
% % NEED TO CHECK THIS PART
% %
% if(m < n)
% % disp('fast svd')
% %
% % Apply A' to a random matrix, obtaining H.
% %
%   if(isreal(A))
%     H = (randn(l,m)*A)';
%   end
%   if(~isreal(A))
%     H = (( randn(l,m) + 1i*randn(l,m) )*A)';
%   end
% %
% % Initialize F to its final size and fill its leftmost block with H.
% %
%   F = zeros(n,(its+1)*l);
%   F(1:n, 1:l) = H;
% %
% % Apply A'*A to H a total of its times,
% % augmenting F with the new H each time.
% %
%   for it = 1:its
%     H = A*H;
%     H = (H'*A)';
%     F(1:n, (1+it*l):((it+1)*l)) = H;
%   end
%   clear H;
% %
% % Form a matrix Q whose columns constitute an orthonormal basis
% % for the columns of F.
% %
%   [Q,~] = qr(F,0);
%   
%   Q = Q(:,1:k);
%   R = Q'*A;
% 
%   clear F;
% 
% end

function [Q,R] = svd_fft(A,l)
% [U,S,V] = svd_fft(A,k,l)
%PCA  Low-rank approximation in SVD form.
%
%
%   [U,S,V] = PCA(A)  constructs a nearly optimal rank-6 approximation
%             USV' to A, using 2 full iterations of a block Lanczos method
%             of block size 6+2=8, started with an n x 8 random matrix,
%             when A is m x n; the ref. below explains "nearly optimal."
%             The smallest dimension of A must be >= 6 when A is
%             the only input to PCA.
%
%   [U,S,V] = PCA(A,k)  constructs a nearly optimal rank-k approximation
%             USV' to A, started with an n x (k+2) random matrix,
%             when A is m x n; the ref. below explains "nearly optimal."
%             k must be a positive integer <= the smallest dimension of A.
%
%
%   [U,S,V] = PCA(A,k,l)  constructs a nearly optimal rank-k approx.
%             USV' to A, using subsampled fft, started with an n x l random matrix,
%             when A is m x n; the ref. below explains "nearly optimal."
%             k must be a positive integer <= the smallest dimension of A,
%             its must be a nonnegative integer,
%             and l must be a positive integer >= k.
%
%
%   The low-rank approximation USV' is in the form of an SVD in the sense
%   that the columns of U are orthonormal, as are the columns of V,
%   the entries of S are all nonnegative, and the only nonzero entries
%   of S appear in non-increasing order on its diagonal.
%   U is m x k, V is n x k, and S is k x k, when A is m x n.
%
%   Increasing its or l improves the accuracy of the approximation USV'
%   to A; the ref. below describes how the accuracy depends on its and l.
%
%
%   Note: PCA invokes RAND. To obtain repeatable results,
%         invoke RAND('seed',j) with a fixed integer j before invoking PCA.
%
%   Note: PCA currently requires the user to center and normalize the rows
%         or columns of the input matrix A before invoking PCA (if such
%         is desired).
%
%   Note: The user may ascertain the accuracy of the approximation USV'
%         to A by invoking DIFFSNORM(A,U,S,V).
%
%
%   inputs (the first is required):
%   A -- matrix being approximated
%   k -- rank of the approximation being constructed;
%        k must be a positive integer <= the smallest dimension of A,
%        and defaults to 6
%   its -- number of full iterations of a block Lanczos method to conduct;
%          its must be a nonnegative integer, and defaults to 2
%   l -- block size of the block Lanczos iterations;
%        l must be a positive integer >= k, and defaults to k+2
%
%   outputs (all three are required):
%   U -- m x k matrix in the rank-k approximation USV' to A,
%        where A is m x n; the columns of U are orthonormal
%   S -- k x k matrix in the rank-k approximation USV' to A,
%        where A is m x n; the entries of S are all nonnegative,
%        and its only nonzero entries appear in nonincreasing order
%        on the diagonal
%   V -- n x k matrix in the rank-k approximation USV' to A,
%        where A is m x n; the columns of V are orthonormal
%
%
%   Example:
%      A = rand(1000,2)*rand(2,1000);
%      A = randn(1000,1000);
%     A = A/normest(A);
%      [U,S,V] = pca2(A,2,0);
%      diffsnorm(A,U,S,V)
%
%     This code snippet produces a rank-2 approximation USV' to A such that
%     the columns of U are orthonormal, as are the columns of V, and
%     the entries of S are all nonnegative and are zero off the diagonal.
%     diffsnorm(A,U,S,V) outputs an estimate of the spectral norm
%     of A-USV', which should be close to the machine precision.
%
%
%   Reference:
%   Vladimir Rokhlin, Arthur Szlam, and Mark Tygert,
%   A randomized algorithm for principal component analysis,
%   arXiv:0809.2274v1 [stat.CO], 2008 (available at http://arxiv.org).
%
%
%   See also PCACOV, PRINCOMP, SVDS.
%

%   Copyright 2009 Mark Tygert.

%
% Retrieve the dimensions of A.
%
[m,n] = size(A);

if(m >= n)
%
% Apply A to a random matrix, obtaining H.
%
  Yaux = ADmult(A);
  Yaux = (fft(Yaux'))';
  [m,n] = size(Yaux);
  clist = randi(n,1,l);
  Y = zeros(m,l);
  for ind = 1:l
      Y(:,ind) = Yaux(:,clist(ind));
  end
  [Q,R] = qr(Y,0);

  clear Y;
end

if(m < n)
%
% Apply A' to a random matrix, obtaining H.
%
  Yaux = ADmult(A');
  Yaux = (fft(Yaux'))';
  [m,n] = size(Yaux);
  clist = randi(n,1,l);
  Y = zeros(m,l);
  for ind = 1:l
      Y(:,ind) = Yaux(:,clist(ind));
  end
  [Q,R] = qr(Y,0);
  clear Y;
end

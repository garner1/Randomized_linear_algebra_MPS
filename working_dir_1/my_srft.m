function C = my_srft(A,l)
%
% Given A of size n, k<=n, and l<=n,
% returns  C = A*S and W = S'*A*S, where S is an n*l SRFT matrix, 
% that is, S = sqrt(n/l)*D F^t R^t, where D is a diagonal matrix of random 
% signs, F is the FFT matrix, and R restricts from n coordinates to l
%

n = size(A,1);
colindices = randperm(n);
colindices = colindices(1:l);

d = 1-2*(rand(1,n) > .5);

%
% fast A*D
%
Y = fft(bsxfun(@times,A,d)'); 
C = Y(colindices, :)';



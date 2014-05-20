function out = ADmult(A)

[m,n] = size(A);
out = zeros(m,n);
list = exp(1i*2*pi*rand(n,1));
for ind = 1:n
    out(:,ind) = list(ind).*A(:,ind);
end


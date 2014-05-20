function out = applyfft(M)

% [m,n] = size(M);
% out = zeros(m,n);
% for ind = 1:m
%     out(ind,:) = fft(M(ind,:));
% end

out = fft(M');
out = out';

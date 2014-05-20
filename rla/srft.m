% PREPARE A SUB-SAMPLED RANDOM FOURIER TRANSFORM
% REF: FINDING STRUCTURE WITH RANDOMNESS
function [D,F,R] = srft(n,l)

D = sparse(1:n,1:n,exp(1i*2*pi*rand(n,1)),n,n,n)/sqrt(l);

F = dftmtx(n);

R = sparse(randi(n,1,l),1:l,ones(1,l),n,l,l);


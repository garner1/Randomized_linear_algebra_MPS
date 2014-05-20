function [mps] = createrandommps(N,D,d)

mps = cell(1,N); 
%
%to be used only with OBC
%
% mps{1}=randn(1,D,d)/sqrt(D); 
% mps{N}=randn(D,1,d)/sqrt(D); 
for i = 1:N
%     for i = 2:N-1 % to be used with OBC
    mps{i}=randn(D,D,d)/sqrt(D); 
end
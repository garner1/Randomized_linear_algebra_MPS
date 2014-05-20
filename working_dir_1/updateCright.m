function [Cright]=updateCright(Cright,B,X,A)

if isempty(X), X=reshape(eye(size(B,3)),[1,1,2,2]); end

%A[dl,dr,s]Cright[i1,i2,i3]
Cright=contracttensors(A,3,2,Cright,3,3); 
Cright=contracttensors(X,4,[2,4],Cright,4,[4,2]); 
Cright=contracttensors(conj(B),3,[2,3],Cright,4,[4,2]);
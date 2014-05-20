function [Cleft]=updateCleft(Cleft,B,X,A)

if isempty(X), X=reshape(eye(size(B,3)),[1,1,2,2]); end

Cleft=contracttensors(A,3,1,Cleft,3,3); 
Cleft=contracttensors(X,4,[1,4],Cleft,4,[4,2]); 
Cleft=contracttensors(conj(B),3,[1,3],Cleft,4,[4,2]);

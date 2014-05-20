% ************************ one-site optimization **************************

function [B,K]=reduceD2_onesite(A,X,Cleft,Cright)

Cleft=contracttensors(Cleft,3,3,A,3,1); 
Cleft=contracttensors(Cleft,4,[2,4],X,4,[1,4]);

B=contracttensors(Cleft,4,[3,2],Cright,3,[2,3]); 
B=permute(B,[1,3,2]);

b=reshape(B,[prod(size(B)),1]); 
K=-b'*b;
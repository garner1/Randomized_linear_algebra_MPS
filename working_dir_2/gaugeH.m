function [Uh,Vh] = gaugeH(Uh,Vh,X,Y)

if isempty(X), X=eye(size(Uh,1));end
if isempty(Y), Y=eye(size(Uh,1));end

Uh = contracttensors(Uh,4,3,X,2,2);
Uh = contracttensors(conj(X),2,2,Uh,4,1);
Uh = permute(Uh,[1,2,4,3]);

Vh = contracttensors(Vh,4,2,Y.',2,1);
Vh = contracttensors(Vh,4,3,Y',2,1);
Vh = permute(Vh,[1,3,2,4]);



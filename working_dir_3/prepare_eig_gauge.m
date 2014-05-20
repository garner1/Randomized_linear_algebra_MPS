function [Heff,Neff,X,Y,V1,U1,V2,U2] = prepare_eig_gauge(V1,mpo1,U1,V2,U2,D,d)

U2 = permute(U2,[1 3 2 4]);
X = sqrtm(pinv(U2));
V2 = permute(V2,[2 4 1 3]);
Y = sqrtm(pinv(V2));
% V2 = Y*V2*Y';
% U2 = X*U2*X';
% fact = V2(1,1)*U2(1,1);
%
%gauge N
%
U2 = contracttensors(X,2,2,U2,2,1);
U2 = contracttensors(U2,2,2,conj(X),2,2);
V2 = contracttensors(Y.',2,1,V2,2,1);
V2 = contracttensors(V2,2,2,Y',2,1);
Neff = contracttensors(V2,3,3,eye(d),3,3);
Neff = contracttensors(Neff,5,5,U2,3,3);
Neff = permute(Neff,[1 3 5 2 4 6]);
Neff= reshape(Neff,[D*d*D,D*d*D]);

V2 = permute(V2,[4 1 3 2]);
U2 = permute(U2,[1 3 2 4]);
%
%gauge the hamiltonian
%
U1 = contracttensors(X,2,2,U1,4,1);
U1 = permute(contracttensors(U1,4,3,conj(X),2,2),[1 2 4 3]);
V1 = contracttensors(Y.',2,1,V1,4,2);
V1 = permute(contracttensors(V1,4,4,Y',2,1),[2 1 3 4]);

Heff = contracttensors(V1,4,3,mpo1,4,1);
Heff = contracttensors(Heff,6,[1,4],U1,4,[4,2]);
Heff = permute(Heff,[1 3 5 2 4 6]);
Heff= reshape(Heff,[D*d*D,D*d*D]);


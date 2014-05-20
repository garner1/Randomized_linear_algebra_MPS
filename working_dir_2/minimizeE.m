function [mpsA,evalue,nn,U] = minimizeE(mpox,Vh,Uh,Un,Vn,X,Y,v0)
%
%Construct the effective Hamiltonian
%
if isempty(X), X=eye(size(Uh,1));end
if isempty(Y), Y=eye(size(Uh,1));end

D = size(Vh,2);
d = size(mpox,3);
Heff = contracttensors(Vh,4,3,mpox,4,1);
Heff = contracttensors(Heff,6,[1,4],Uh,4,[4,2]);
Heff = permute(Heff,[1,5,3,2,6,4]);
Heff = reshape(Heff,[D^2*d,D^2*d]);
opts.v0 = v0;
[aux,evalue] = eigs(Heff,1,'sr',opts);
aux = reshape(aux,[D,D,d]);
%
%regauge
%
aux = contracttensors(Y.',2,2,aux,3,1);
aux = contracttensors(aux,3,2,X,2,1);
mpsA = permute(aux,[1,3,2]);
%
%norm
%
nn = norma(mpsA,Un,Vn);
%
%canonical form
%
[mpsA,U] = prepare_onesite(mpsA,'lr');
%
%energy normalized
%
evalue = real(evalue)/real(nn);

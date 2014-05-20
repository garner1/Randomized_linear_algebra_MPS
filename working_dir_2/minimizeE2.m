function [mpsA,evalue,nn,U] = minimizeE2(mpoA,Vh,Uh,Un,Vn,X,Y,v0)
%
%Construct the effective Hamiltonian
%
if isempty(X), X=eye(size(Uh,1));end
if isempty(Y), Y=eye(size(Uh,1));end

D = size(Vh,2);
d = size(mpoA,3);
Heff = contracttensors(Vh,4,3,mpoA,4,1);
Heff = contracttensors(Heff,6,[1,4],Uh,4,[4,2]);
Heff = permute(Heff,[1,5,3,2,6,4]);
Heff = reshape(Heff,[D^2*d,D^2*d]);

Neff = contracttensors(Vn,4,3,reshape(eye(2),[1,1,2,2]),4,1);
Neff = contracttensors(Neff,6,[1,4],Un,4,[4,2]);
Neff = permute(Neff,[1,5,3,2,6,4]);
Neff = reshape(Neff,[D^2*d,D^2*d]);

opts.v0 = v0;
[aux,evalue] = eigs(Heff,Neff,1,'sr',opts);

mpsA = reshape(aux,[D,D,d]);
%
%regauge
%main
mpsA = contracttensors(Y.',2,2,mpsA,3,1);
mpsA = contracttensors(mpsA,3,2,X,2,1);
mpsA = permute(mpsA,[1,3,2]);
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
evalue = real(evalue);

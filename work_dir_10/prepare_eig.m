function [Heff,Neff] = prepare_eig(V1,mpo1,U1,V2,mpo2,U2,D,d)

Heff = contracttensors(V1,4,3,mpo1,4,1);
Heff = contracttensors(Heff,6,[1,4],U1,4,[4,2]);
Heff = permute(Heff,[1 3 5 2 4 6]);
Heff= reshape(Heff,[D*d*D,D*d*D]);

Neff = contracttensors(V2,4,3,mpo2,4,1);
Neff = contracttensors(Neff,6,[1,4],U2,4,[4,2]);
Neff = permute(Neff,[1 3 5 2 4 6]);
Neff= reshape(Neff,[D*d*D,D*d*D]);
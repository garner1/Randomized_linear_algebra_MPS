function [mps,evalue,U] = random_eig(mpox,Vh,Uh,X,Y,peig)
% 
%perform fast eig of Heff
%
if isempty(X), X=eye(size(Uh,1));end
if isempty(Y), Y=eye(size(Uh,1));end

D = size(Vh,2);
d = size(mpox,3);
Yeig = zeros(D^2*d,peig);
for col = 1:peig
    p0 = randn(D,D,d)+1i*randn(D,D,d); 
    aux = contracttensors(Vh,4,3,mpox,4,1);
    aux = contracttensors(aux,6,[3,6],p0,3,[1,3]);
    aux = contracttensors(aux,5,[5,3,1],Uh,4,[3,2,4]);
    aux = permute(aux,[1,3,2]);
    Yeig(:,col) = reshape(aux,[D*D*d,1]);
end
[qeig,~] = qr(Yeig,0);
Beig = zeros(peig,peig);
for col = 1:peig
    for row = 1:peig
        p0 = reshape(qeig(:,col),[D,D,d]);
        aux = contracttensors(Vh,4,3,mpox,4,1);
        aux = contracttensors(aux,6,[3,6],p0,3,[1,3]);
        aux = contracttensors(aux,5,[5,3,1],Uh,4,[3,2,4]);
        aux = permute(aux,[1,3,2]);
        Beig(row,col) = contracttensors(conj(reshape(qeig(:,row),[D,D,d])),3,[1,2,3],aux,3,[1,2,3]);
    end
end
[Veig,Deig] = eig(Beig);
aux = reshape(qeig*Veig(:,1),[D,D,d]);
aux = contracttensors(Y.',2,2,aux,3,1);
aux = contracttensors(aux,3,2,X,2,1);
aux = permute(aux,[1,3,2]);
[mps,U] = prepare_onesite(aux,'lr');
evalue = real(Deig(1,1));

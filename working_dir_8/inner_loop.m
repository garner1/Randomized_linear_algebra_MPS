function [mpsma,mpsmb,ee] = inner_loop(V1,mpo1a,mpo1b,U1,V2,mpo2,U2,mpsma,mpsmb,e0,m)
%
%inner loop
%
v = cell(m+1,1); 
[D,~,d] = size(mpsma);

v{1} = permute(contracttensors(mpsma,3,2,mpsmb,3,1),[1,3,2,4]);
v{1} = reshape(v{1},[D,D,d*d]);

mpo1 = contracttensors(mpo1a,4,2,mpo1b,4,1);
mpo1 = permute(mpo1,[1,4,2,5,3,6]);
mpo1 = reshape(mpo1,[5,5,d*d,d*d]);

mpo2 = contracttensors(mpo2,4,2,mpo2,4,1);
mpo2 = permute(mpo2,[1,4,2,5,3,6]);
mpo2 = reshape(mpo2,[1,1,d*d,d*d]);

vv = zeros(D*D*d*d,m+1);
vv(:,1) = reshape(v{1},[D*D*d*d,1]);
for ii = 1:m
    v1 = matvectn(V1,U1,mpo1,v{ii});
    v2 = matvectn(V2,U2,mpo2,v{ii});
    v{ii+1} = v1-e0*v2;
    vv(:,ii+1) = reshape(v{ii+1},[D*D*d*d,1]);
end
%
%find basis
%
[zz,~] = qr(vv,0);
%
%form the small pencil
%
mm = size(zz,2);
Am = zeros(mm);
Bm = zeros(mm);
for jj = 1:mm
    z = reshape(zz(:,jj),[D,D,d*d]);
    v1 = matvectn(V1,U1,mpo1,z);
    v2 = matvectn(V2,U2,mpo2,z);
    zcol = v1-e0*v2;
    for ii = 1:mm
        zrow = reshape(zz(:,ii)',[D,D,d*d]);
        Am(ii,jj) = contracttensors(zrow,3,[1,2,3],zcol,3,[1,2,3]);
        Bm(ii,jj) = contracttensors(zrow,3,[1,2,3],v2,3,[1,2,3]);
    end
end
[nu,mu] = eigs(Am,Bm,1,'sr');
ee = e0+mu;
mpsm = reshape(zz*nu,[D,D,d,d]);

mpsm = reshape(permute(mpsm,[1,3,2,4]),[D*d,D*d]);
[U,S,V] = svd2(mpsm);
U = U(:,1:D);
V = S(1:D,1:D)*V(1:D,:);
mpsma = permute(reshape(U,[D,d,D]),[1 3 2]);
mpsmb = permute(reshape(V,[D,d,D]),[1 3 2]);


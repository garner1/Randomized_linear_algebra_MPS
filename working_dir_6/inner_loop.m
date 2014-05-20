function [mpsm,ee] = inner_loop(V1,mpo1,U1,V2,mpo2,U2,mpsm,e0,m)
%
%inner loop
%
v = cell(m+1,1); 
[D,~,d] = size(mpsm);

v{1} = mpsm;
vv = zeros(D*D*d,m+1);
vv(:,1) = reshape(v{1},[D*D*d,1]);
for ii = 1:m
    v1 = matvectn(V1,U1,mpo1,v{ii});
    v2 = matvectn(V2,U2,mpo2,v{ii});
    v{ii+1} = v1-e0*v2;
    vv(:,ii+1) = reshape(v{ii+1},[D*D*d,1]);
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
    z = reshape(zz(:,jj),[D,D,d]);
    v1 = matvectn(V1,U1,mpo1,z);
    v2 = matvectn(V2,U2,mpo2,z);
    zcol = v1-e0*v2;
    for ii = 1:mm
        zrow = reshape(zz(:,ii)',[D,D,d]);
        Am(ii,jj) = contracttensors(zrow,3,[1,2,3],zcol,3,[1,2,3]);
        Bm(ii,jj) = contracttensors(zrow,3,[1,2,3],v2,3,[1,2,3]);
    end
end
[nu,mu] = eigs(Am,Bm,1,'sr');
ee = e0+mu;
mpsm = reshape(zz*nu,[D,D,d]);
        

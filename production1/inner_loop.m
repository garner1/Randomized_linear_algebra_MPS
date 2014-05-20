function [mpsm,ee] = inner_loop(V1,mpo1,U1,V2,mpo2,U2,mpsm,e0,m)
%
%inner loop
%
[D,~,d] = size(mpsm);

vv = zeros(D*D*d,m+1);
z1 = zeros(D,D,d,m);
vv(:,1) = reshape(mpsm,[D*D*d,1]);
for ii = 1:m
    z1(:,:,:,ii) = reshape(vv(:,ii),[D,D,d]);
    vv(:,ii+1) = reshape(matvectn(V1,U1,mpo1,z1(:,:,:,ii)),[D*D*d,1])...
        -e0*reshape(matvectn(V2,U2,mpo2,z1(:,:,:,ii)),[D*D*d,1]);
end
%
%find basis
%
[zz,~] = qr(vv,0);
%
%form the small pencil
%
mm = size(zz,2);
vv2 = zeros(D*D*d,mm);v2 = vv2;
z2 = zeros(D,D,d,m);
for jj = 1:mm
    z2(:,:,:,ii) = reshape(zz(:,jj),[D,D,d]);
    v2(:,jj) = reshape(matvectn(V2,U2,mpo2,z2(:,:,:,ii)),[D*D*d,1]);
    vv2(:,jj) = reshape(matvectn(V1,U1,mpo1,z2(:,:,:,ii)),[D*D*d,1])-e0*v2(:,jj);
end
Am = zz'*vv2;
Bm = zz'*v2;
[nu,mu] = eigs(Am,Bm,1,'sr');
ee = real(e0+mu);
mpsm = reshape(zz*nu,[D,D,d]);


clear all;
N = 2^6;
D = 2^1;
d = 2;

jx = 1;
jy = 0;
jz = 0;
hz = 1; h = hz;
%
%prepare mpo
%
mpo = mpo_hh(N,jx,jy,jz,h);
mpoid = mpo_id(N,d); 
%
%prepare mps
%
mps = randmps_norm(N,D);

ppp = 5;
iter = 0;
% l = 2+ppp;
for site = 1:1
    site=15;
    %
    %Perform random svd for the normalization
    %
    r = ppp+3;
    [U,S,V,listn,j] = adaptive_svd(mps,mpoid,site,ppp,iter,r,1e-6);
    [p,pp,error] = exact(mpoid,mps,U,V,site,D,ppp);
    [u,s,v]=svd(p);

    aux = reshape(permute(contracttensors(conj(mps{site}),3,3,mps{site},3,3),[1,3,2,4]),[D*D,D*D]);
    P = reshape(U,[D*1*D,ppp])*reshape(V,[ppp,D*1*D]);
    PP = reshape(permute(contracttensors(U,4,[2,4],V,4,[3,1]),[1,3,2,4]),[D*D,D*D]);

    disp(error);
    disp(abs(trace((P-p)*aux)/trace(p*aux)));
    disp('---');
end
%main
%test3 exact

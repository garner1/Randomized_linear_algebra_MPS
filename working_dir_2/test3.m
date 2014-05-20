clear all;
N = 2^7;
D = 2^3;
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
ppp=4;
iter=0;
l=2+ppp;
for site = 1:N
    %
    %Perform random svd for the normalization
    %
    [U,S,V] = random_svd(mps,mpoid,site,ppp,iter,l);
    [p,~,error] = exact(mpoid,mps,U,V,site,D,ppp);
%     [u,s,v]=svd(p);
    disp(error);
%     disp(abs(error-s(ppp+1,ppp+1)));
    disp('---');
%         r = pn+3;
%         [Un,Sn,Vn] = adaptive_svd(mps,mpoid,site,pn,0,r,1e-12);
%         [X,Y] = gaugeN(Un,Vn);
end
%main
%test3 exact

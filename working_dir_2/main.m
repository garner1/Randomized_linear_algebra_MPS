%test3
clear all;
N = 2^5;
D = 2^2;
d = 2;
sweeps = 1;

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

elist1 = [];elist2 = [];nlist = [];slist = [];
count = 0; 
iter = 2;
pn = 1;ln = pn+2;
ph = 4;lh = 2*ph;
while count<=sweeps,
    count = count+1;
    for site = 1:N
%         p = p+2; l = p+2;
%         if count == 1, p = 1; l = 2*p;end
%         if count == 2, p = 5; l = p*2;end
%         if count == 3, p = 5; l = p*2;end
%         if count == 4, p = 5; l = p*2;end
        %
        %Perform random svd for the normalization
        %
        [Un,Sn,Vn] = random_svd(mps,mpoid,site,pn,iter,ln);
%         r = pn+3;
%         [Un,Sn,Vn] = adaptive_svd(mps,mpoid,site,pn,0,r,1e-12);
%         [X,Y] = gaugeN(Un,Vn);
        %
        %Perform random svd for the effective Hamiltonian
        %
        [Uh,Sh,Vh] = random_svd(mps,mpo,site,ph,iter,lh);
%         r = ph+2;
%         [Uh,Sh,Vh] = adaptive_svd(mps,mpo,site,ph,0,r,1e-10);
%         [Uh,Vh] = gaugeH(Uh,Vh,[],[]);

        [mps{site},evalue2,nn2,U] = minimizeE2(mpo{site},Vh,Uh,Un,Vn,[],[],reshape(mps{site},[D^2*d,1]));
        if site==N,
            mps{site+1} = contracttensors(U,2,2,mps{1},3,1);
        else
            mps{site+1} = contracttensors(U,2,2,mps{site+1},3,1);
        end
        
        
        disp([count,site]);
        disp(real(evalue2)/N);
%         disp(real(nn2));
        disp('---')
        elist2 = [elist2,real(evalue2)/N];
        plot(1:N*(count-1)+site,elist2);drawnow;
    end
end
%         [mps{site},evalue2,nn2,U] = minimizeE(mpo{site},Vh,Uh,Un,Vn,[],[],reshape(mps{site},[D^2*d,1]));
%         disp([count,site,real(evalue1)/N,real(nn1),real(evalue2)/N,real(nn2)]);
%         elist1 = [elist1,real(evalue1)/N];
%         plot(1:N*(count-1)+site,elist1,1:N*(count-1)+site,elist2);drawnow;

% clear all;
N = 2*2^9;
D = 2^5;
d = 2;
sweeps = 2;

pn = 1; ln = 2+pn;
ph = 2; lh = 2+ph;
iter = 0;

jx = 1;
jy = 0;
jz = 0;
hz = 1; h = hz;

Vn = cell(1,N);
Vh = cell(1,N);
mpo = mpo_hh(N,jx,jy,jz,h);
mpoid = mpo_id(N,d); 
mps = randmps_norm(N,D);

elist = [];
count = 0;
while count < sweeps, 
    count = count+1;
    %
    %solve first cluster
    %
    [Un_env,Sn,Vn_env] = randomEnv_svd(mps(1,N/2+1:N),mpoid(1,N/2+1:N),pn,iter,ln);
    [Uh_env,Sh,Vh_env] = randomEnv_svd(mps(1,N/2+1:N),mpo(1,N/2+1:N),ph,iter,lh);
    Vh{1} = Vh_env;
    Vn{1} = Vn_env;
    Uh = init_renv(mps,mpo,Uh_env,'two');
    Un = init_renv(mps,mpoid,Un_env,'two');
    for site = 1:N/2
        [Heff,Neff,X,Y,Vh{site},Uh{site},Vn{site},Un{site}] = prepare_eig_gauge(Vh{site},mpo{site},Uh{site},Vn{site},Un{site},D,d);
        opts.v0 = reshape(permute(mps{site},[1 3 2]),[D*d*D,1]);
        [aux,evalue] = eigs(Heff,1,'sr',opts);
        mps{site} = permute(reshape(aux,[D,d,D]),[1,3,2]);
        
        [mps{site},U] = prepare_onesite(mps{site},'lr');
        if site==N/2,
            mps{1+N/2} = contracttensors(U,2,2,mps{1+N/2},3,1);
        else
            mps{site+1} = contracttensors(U,2,2,mps{site+1},3,1);
            Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
            Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
        end
        elist = [elist,real(evalue)/N];
        plot(elist,'-o');drawnow;
        disp([count,site,real(evalue)/N])
% DD=real(eig(Neff));
% disp(norm(Neff-eye(size(Neff))));
    end
    %
    %solve second cluster
    %
    [Un_env,Sn,Vn_env] = randomEnv_svd(mps(1,1:N/2),mpoid(1,1:N/2),pn,iter,ln);
    [Uh_env,Sh,Vh_env] = randomEnv_svd(mps(1,1:N/2),mpo(1,1:N/2),ph,iter,lh);
    Vh{N/2+1} = Vh_env;
    Vn{N/2+1} = Vn_env;
    Uh = init_renv(mps,mpo,Uh_env,'one');
    Un = init_renv(mps,mpoid,Un_env,'one');
    for site = N/2+1:N
        [Heff,Neff,X,Y,Vh{site},Uh{site},Vn{site},Un{site}] = prepare_eig_gauge(Vh{site},mpo{site},Uh{site},Vn{site},Un{site},D,d);
        opts.v0 = reshape(permute(mps{site},[1 3 2]),[D*d*D,1]);
        [aux,evalue] = eigs(Heff,1,'sr',opts);
        mps{site} = permute(reshape(aux,[D,d,D]),[1,3,2]);
        
        [mps{site},U] = prepare_onesite(mps{site},'lr');
        if site==N,
            mps{1} = contracttensors(U,2,2,mps{1},3,1);
        else
            mps{site+1} = contracttensors(U,2,2,mps{site+1},3,1);
            Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
            Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
        end
        elist = [elist,real(evalue)/N];
        plot(elist,'-o');drawnow;
        disp([count,site,real(evalue)/N])
% DD=real(eig(Neff));
% disp(norm(Neff-eye(size(Neff))));
    end
end


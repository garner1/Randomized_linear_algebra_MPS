clear all;
N = 160;
D = 2^5;
size_bulk= 40;
d = 2;
sweeps = round(4*N/size_bulk);
pn = 2^1; ln = 2*pn; %used in rand approx
ph = 2^2; lh = 2*ph; %used in rand approx
iter = 2; %used in rand approx
m_in = 10; %used in the inner loop of the solver
m_out = 5; %used in the outer loop of the solver
precision = 1*1e-4; %used in the outer loop of the solver
        
phi = 0.0;
J = 0.25;
jx = 1.0*J;
jy = 1.0*J;
jz = 0.0*J;
hz = 0; h = hz;

Ln = cell(1,N);
Lh = cell(1,N);
mpo = mpo_hh(N,jx,jy,jz,h,phi);
mpoid = mpo_id(N,d); 
mps = randmps_norm(N,D);

elist = [];
elist2 = [];
count = 0;
site_in = 1;
site_den = 0;
while count < sweeps, 
    count = count+1;
    site_fin = site_in+size_bulk-1;
    if site_fin > N, site_fin = mod(site_fin,N);end
    if site_in > N, site_in = mod(site_in,N);end
    if site_fin>site_in,
tic
        [Rn_env,Ln_env] = randomEnv(mps(1,[site_fin+1:N,1:site_in-1]),mpoid(1,[site_fin+1:N,1:site_in-1]),pn,iter,ln);
toc
tic
        [Rh_env,Lh_env] = randomEnv(mps(1,[site_fin+1:N,1:site_in-1]),mpo(1,[site_fin+1:N,1:site_in-1]),ph,iter,lh);
toc
        Lh{site_in} = Lh_env;
        Ln{site_in} = Ln_env;
        Rh = init_renv(mps(1,site_in:site_fin),mpo(1,site_in:site_fin),Rh_env);
        Rn = init_renv(mps(1,site_in:site_fin),mpoid(1,site_in:site_fin),Rn_env);

        sitex = 0;
        for site = site_in:site_fin
            site_den = site_den+1;
            sitex = sitex + 1;
            [aux,evalue,list] = outer_loop(Lh{site},mpo{site},Rh{sitex},Ln{site},mpoid{site},Rn{sitex},mps{site},m_in,m_out,precision,N);            
%             value = expvalue(mpo{site},aux,Lh{site},Rh{sitex},Ln{site},Rn{sitex});
            elist = [elist,evalue/N];
            mps{site} = aux;
            [mps{site},U] = prepare_onesite(mps{site},'lr');
            if site == N,
                mps{1} = contracttensors(U,2,2,mps{1},3,1);
                Lh{site+1} = update_lenv(mps{site},mpo{site},Lh{site});
                Ln{site+1} = update_lenv(mps{site},mpoid{site},Ln{site});
            elseif site==site_fin,
                mps{1+site_fin} = contracttensors(U,2,2,mps{1+site_fin},3,1);
            else
                mps{site+1} = contracttensors(U,2,2,mps{site+1},3,1);
                Lh{site+1} = update_lenv(mps{site},mpo{site},Lh{site});
                Ln{site+1} = update_lenv(mps{site},mpoid{site},Ln{site});
            end
%             elist2 = [elist2,real(value)/N];
            plot(elist);drawnow;
            disp([count,site,elist(end)])
        end
        
    else
tic
        [Rn_env,Ln_env] = randomEnv2(mps(1,site_fin+1:site_in-1),mpoid(1,site_fin+1:site_in-1),pn,iter,ln);
toc
tic
        [Rh_env,Lh_env] = randomEnv2(mps(1,site_fin+1:site_in-1),mpo(1,site_fin+1:site_in-1),ph,iter,lh);
toc        
        Lh{site_in} = Lh_env;
        Ln{site_in} = Ln_env;
        Rh = init_renv(mps(1,[site_in:N,1:site_fin]),mpo(1,[site_in:N,1:site_fin]),Rh_env);
        Rn = init_renv(mps(1,[site_in:N,1:site_fin]),mpoid(1,[site_in:N,1:site_fin]),Rn_env);
        
        sitex = 0;
        for site = [site_in:N,1:site_fin]
            site_den = site_den+1;
            sitex = sitex+1;
            [aux,evalue,list] = outer_loop(Lh{site},mpo{site},Rh{sitex},Ln{site},mpoid{site},Rn{sitex},mps{site},m_in,m_out,precision,N);            
%             value = expvalue(mpo{site},aux,Lh{site},Rh{sitex},Ln{site},Rn{sitex});
            elist = [elist,evalue/N];
            mps{site} = aux;
            [mps{site},U] = prepare_onesite(mps{site},'lr');
            if site==site_fin,
                mps{1+site_fin} = contracttensors(U,2,2,mps{1+site_fin},3,1);
            elseif site==N,
                mps{1} = contracttensors(U,2,2,mps{1},3,1);
                Lh{1} = update_lenv(mps{site},mpo{site},Lh{site});
                Ln{1} = update_lenv(mps{site},mpoid{site},Ln{site});
            else
                mps{site+1} = contracttensors(U,2,2,mps{site+1},3,1);
                Lh{site+1} = update_lenv(mps{site},mpo{site},Lh{site});
                Ln{site+1} = update_lenv(mps{site},mpoid{site},Ln{site});
            end
%             elist2 = [elist2,evalue/N];
            plot(elist);drawnow;
            disp([count,site,elist(end)])
        end
    end          
    site_in = site_in+round(N/size_bulk-1)*size_bulk;
end


clear all;
N = 3*2^8;
D = 2^6;
d = 2;
sweeps = 1;
m_in = 3; 
m_out = 3;

jx = 0.25;
jy = 0.25;
jz = 0.25;
hz = 0; h = hz;

Vn = cell(1,N);
Vh = cell(1,N);
mpo = mpo_hh(N,jx,jy,jz,h);
mpoid = mpo_id(N,d); 
mps = randmps_norm(N,D);

elist = [];
count = 0;
while count < sweeps, 
    count = count+1;
    if count==1,
        pn = 2^0; ln = 2+pn;
        ph = 2^1; lh = 2+ph;
        iter = 0;
    end
    if count>=2,
        pn = 2^0; ln = 2+pn;
        ph = 2^1; lh = 2+ph;
        iter = 0;
    end
    %
    %solve first cluster
    %
    [Un_env2,Vn_env2] = randomEnv(mps(1,N/3+1:2*N/3),mpoid(1,N/3+1:2*N/3),pn,iter,ln);
    [Uh_env2,Vh_env2] = randomEnv(mps(1,N/3+1:2*N/3),mpo(1,N/3+1:2*N/3),ph,iter,lh);
    [Un_env3,Vn_env3] = randomEnv(mps(1,2*N/3+1:N),mpoid(1,2*N/3+1:N),pn,iter,ln);
    [Uh_env3,Vh_env3] = randomEnv(mps(1,2*N/3+1:N),mpo(1,2*N/3+1:N),ph,iter,lh);
    Vh{1} = Vh_env3;
    Vn{1} = Vn_env3;
    aux = contracttensors(Vn_env2,4,[2,3,4],Un_env3,4,[1 2 3]);
    Un_env2 = contracttensors(Un_env2,4,4,aux,2,1);
    aux = contracttensors(Vh_env2,4,[2,3,4],Uh_env3,4,[1 2 3]);
    Uh_env2 = contracttensors(Uh_env2,4,4,aux,2,1);
    Uh = init_renv(mps,mpo,Uh_env2,'one');
    Un = init_renv(mps,mpoid,Un_env2,'one');
    for site = 1:N/3
        [aux,evalue] = outer_loop(Vh{site},mpo{site},Uh{site},Vn{site},mpoid{site},Un{site},mps{site},m_in,m_out);            
        aux = permute(reshape(aux,[D,d,D]),[1,3,2]);
        value = expvalue(mpo{site},aux,Vh{site},Uh{site},Vn{site},Un{site});
        if ~isempty(elist) && value/N > elist(end),
            elist = [elist,elist(end)]; 
            Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
            Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
        else
            elist = [elist,value/N];
            mps{site} = aux;
            [mps{site},U] = prepare_onesite(mps{site},'lr');
            if site==N/3,
                mps{1+N/3} = contracttensors(U,2,2,mps{1+N/3},3,1);
            else
                mps{site+1} = contracttensors(U,2,2,mps{site+1},3,1);
                Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
                Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
            end
        end
        plot(elist);drawnow;
        disp([count,site,elist(end)])
    end
    %
    %solve second cluster
    %
    [Un_env1,Vn_env1] = randomEnv(mps(1,1:N/3),mpoid(1,1:N/3),pn,iter,ln);
    [Uh_env1,Vh_env1] = randomEnv(mps(1,1:N/3),mpo(1,1:N/3),ph,iter,lh);
    [Un_env3,Vn_env3] = randomEnv(mps(1,2*N/3+1:N),mpoid(1,2*N/3+1:N),pn,iter,ln);
    [Uh_env3,Vh_env3] = randomEnv(mps(1,2*N/3+1:N),mpo(1,2*N/3+1:N),ph,iter,lh);
    Vh{N/3+1} = Vh_env1;
    Vn{N/3+1} = Vn_env1;
    aux = contracttensors(Vn_env3,4,[2,3,4],Un_env1,4,[1 2 3]);
    Un_env3 = contracttensors(Un_env3,4,4,aux,2,1);
    aux = contracttensors(Vh_env3,4,[2,3,4],Uh_env1,4,[1 2 3]);
    Uh_env3 = contracttensors(Uh_env3,4,4,aux,2,1);
    Uh = init_renv(mps,mpo,Uh_env3,'two');
    Un = init_renv(mps,mpoid,Un_env3,'two');
    for site = N/3+1:2*N/3
        [aux,evalue] = outer_loop(Vh{site},mpo{site},Uh{site},Vn{site},mpoid{site},Un{site},mps{site},m_in,m_out);            
        aux = permute(reshape(aux,[D,d,D]),[1,3,2]);
        value = expvalue(mpo{site},aux,Vh{site},Uh{site},Vn{site},Un{site});
        if ~isempty(elist) && value/N > elist(end),
            elist = [elist,elist(end)]; 
            Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
            Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
        else
            elist = [elist,value/N];
            mps{site} = aux;
            [mps{site},U] = prepare_onesite(mps{site},'lr');
            if site==2*N/3,
                mps{2*N/3+1} = contracttensors(U,2,2,mps{2*N/3},3,1);
            else
                mps{site+1} = contracttensors(U,2,2,mps{site+1},3,1);
                Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
                Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
            end
        end
        plot(elist);drawnow;
        disp([count,site,elist(end)])
    end
    %
    %solve third cluster
    %
    [Un_env2,Vn_env2] = randomEnv(mps(1,N/3+1:2*N/3),mpoid(1,N/3+1:2*N/3),pn,iter,ln);
    [Uh_env2,Vh_env2] = randomEnv(mps(1,N/3+1:2*N/3),mpo(1,N/3+1:2*N/3),ph,iter,lh);
    [Un_env1,Vn_env1] = randomEnv(mps(1,1:N/3),mpoid(1,1:N/3),pn,iter,ln);
    [Uh_env1,Vh_env1] = randomEnv(mps(1,1:N/3),mpo(1,1:N/3),ph,iter,lh);
    Vh{2*N/3+1} = Vh_env2;
    Vn{2*N/3+1} = Vn_env2;
    aux = contracttensors(Vn_env1,4,[2,3,4],Un_env2,4,[1 2 3]);
    Un_env1 = contracttensors(Un_env1,4,4,aux,2,1);
    aux = contracttensors(Vh_env1,4,[2,3,4],Uh_env2,4,[1 2 3]);
    Uh_env1 = contracttensors(Uh_env1,4,4,aux,2,1);
    Uh = init_renv(mps,mpo,Uh_env1,'three');
    Un = init_renv(mps,mpoid,Un_env1,'three');
    for site = 2*N/3+1:N
        [aux,evalue] = outer_loop(Vh{site},mpo{site},Uh{site},Vn{site},mpoid{site},Un{site},mps{site},m_in,m_out);            
        aux = permute(reshape(aux,[D,d,D]),[1,3,2]);
        value = expvalue(mpo{site},aux,Vh{site},Uh{site},Vn{site},Un{site});
        if ~isempty(elist) && value/N > elist(end),
            elist = [elist,elist(end)]; 
            Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
            Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
        else
            elist = [elist,value/N];
            mps{site} = aux;
            [mps{site},U] = prepare_onesite(mps{site},'lr');
            if site==N,
                mps{1} = contracttensors(U,2,2,mps{1},3,1);
            else
                mps{site+1} = contracttensors(U,2,2,mps{site+1},3,1);
                Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
                Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
            end
        end
        plot(elist);drawnow;
        disp([count,site,elist(end)])
    end
end


clear all;
N = 3*2^5;
D = 2^3;
d = 2;
sweeps = 10;

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
    if count==1,
        pn = 2^2; ln = 2*pn;
        ph = 2^3; lh = 2*ph;
        m_in = 2*ph;
        m_out= 5;
        iter = 0;
    end
    if count>=2,
        pn = 2^2; ln = 2*pn;
        ph = 2^3; lh = 2*ph;
        m_in = 2*ph;
        m_out= 5;
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

    if count==1, 
            num = contracttensors(conj(mps{1}),3,[1,2,3],matvectn(Vh{1},Uh{1},mpo{1},mps{1}),3,[1,2,3]);
            den = contracttensors(conj(mps{1}),3,[1,2,3],matvectn(Vn{1},Un{1},mpoid{1},mps{1}),3,[1,2,3]);
            ee = real(num/den);
            elist = [elist,ee/N];
    end
    for site = 1:N/3-1
        [mps{site},mps{site+1},evalue] = ...
        outer_loop(Vh{site},mpo{site},mpo{site+1},Uh{site+1},Vn{site},mpoid{site},Un{site+1},mps{site},mps{site+1},elist(end),m_in,m_out);
        if site==N/3-1,
            [mps{site+1},U] = prepare_onesite(mps{site+1},'lr');
            mps{site+2} = contracttensors(U,2,2,mps{site+2},3,1);
        else
            Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
            Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
        end
        elist = [elist,real(evalue)/N];
        plot(elist);drawnow;
        disp([count,site,real(evalue)/N])
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
    for site = N/3+1:2*N/3-1
        [mps{site},mps{site+1},evalue] = ...
        outer_loop(Vh{site},mpo{site},mpo{site+1},Uh{site+1},Vn{site},mpoid{site},Un{site+1},mps{site},mps{site+1},elist(end),m_in,m_out);
        if site==2*N/3-1,
            [mps{site+1},U] = prepare_onesite(mps{site+1},'lr');
            mps{site+2} = contracttensors(U,2,2,mps{site+2},3,1);
        else
            Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
            Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
        end
        elist = [elist,real(evalue)/N];
        plot(elist);drawnow;
        disp([count,site,real(evalue)/N])
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
    for site = 2*N/3+1:N-1
        [mps{site},mps{site+1},evalue] = ...
        outer_loop(Vh{site},mpo{site},mpo{site+1},Uh{site+1},Vn{site},mpoid{site},Un{site+1},mps{site},mps{site+1},elist(end),m_in,m_out);
        if site==N-1,
            [mps{site+1},U] = prepare_onesite(mps{site+1},'lr');
            mps{1} = contracttensors(U,2,2,mps{1},3,1);
        else
            Vh{site+1} = update_lenv(mps{site},mpo{site},Vh{site});
            Vn{site+1} = update_lenv(mps{site},mpoid{site},Vn{site});
        end
        elist = [elist,real(evalue)/N];
        plot(elist);drawnow;
        disp([count,site,real(evalue)/N])
    end
end


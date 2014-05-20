clear all;
N = 2*2^5;
D = 2^4;
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
%
%determine the environment
%
    [Un_env,Vn_env] = randomEnv(mps(1,N/2+1:N),mpoid(1,N/2+1:N),pn,iter,ln);
    [Uh_env,Vh_env] = randomEnv(mps(1,N/2+1:N),mpo(1,N/2+1:N),ph,iter,lh);
    Vh{1} = Vh_env;
    Vn{1} = Vn_env;
    Uh = init_renv(mps,mpo,Uh_env,'two');
    Un = init_renv(mps,mpoid,Un_env,'two');
    for site = 1:N/2
%
%note that aux is not an eigenstate in general (because of approx method)
%
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
            if site==N/2,
                mps{1+N/2} = contracttensors(U,2,2,mps{1+N/2},3,1);
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
    [Un_env,Vn_env] = randomEnv(mps(1,1:N/2),mpoid(1,1:N/2),pn,iter,ln);
    [Uh_env,Vh_env] = randomEnv(mps(1,1:N/2),mpo(1,1:N/2),ph,iter,lh);
    Vh{N/2+1} = Vh_env;
    Vn{N/2+1} = Vn_env;
    Uh = init_renv(mps,mpo,Uh_env,'one');
    Un = init_renv(mps,mpoid,Un_env,'one');
    for site = N/2+1:N
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

%         if abs(site-1)<=2||abs(site-N/2)<=2, 
%             [Heff,Neff] = prepare_eig(Vh{site},mpo{site},Uh{site},Vn{site},mpoid{site},Un{site},D,d);
%             opts.v0 = reshape(permute(mps{site},[1 3 2]),[D*d*D,1]);
%             [aux,evalue] = eigs(Heff,Neff,1,'sr',opts);
%             aux = permute(reshape(aux,[D,d,D]),[1,3,2]);
%             value = expvalue(mpo{site},aux,Vh{site},Uh{site},Vn{site},Un{site});
%         else

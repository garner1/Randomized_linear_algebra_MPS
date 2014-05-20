clear all;
N = 4;Dh = 2^1;Dv = (2^1)*Dh;
T = 1;
precision = 1e-2;
mps = generate_bulk(N,Dh);
mpsket = cell(1,N+2);
mpsbra = cell(1,N+2);
for ind = 2:N+1
    mpsket{ind} = mps{ind-1};
    mpsbra{ind} = mps{ind-1};
end
M = 4;
hset=cell(M,N+2);
jx=1; jy=1; jz=1;
sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2); 
for m = 1:M, 
    for j = 1:N+2, 
        if j==1, hset{m,j} = eye(Dv,Dv);
        elseif j==N+2, hset{m,j} = eye(Dv,Dv);
        else hset{m,j} = eye(2); 
        end
    end; 
end
hset{1,2} = jx*sx; 
hset{1,3} = sx; 
hset{2,3} = jx*sx; 
hset{2,4} = sx; 
hset{3,4} = jx*sx; 
hset{3,5} = sx; 
num = []; den = num; E=[];ave=[];rerr=[];
for ind = 1:T
    [A,B,C] = rand_boundary(Dh,Dv,1); 
    mpsket{1} = permute(conj(A),[2 1 3]);
    mpsket{N+2} = A;
    mpsbra{1} = permute(conj(C),[2 1 3]);
    mpsbra{N+2} = C;
    hset{4,N+1} = sx; 
    hset{4,N+2} = reshape(B,[Dv Dv]); 
    hset{4,1} = conj(hset{1,N+2}); 
    hset{4,2} = jx*sx; 
    [E0,ket,bra] = minimizeE(mpsbra,hset,mpsket,precision); 
    E = [E,real(E0)];
    ave = [ave,mean(E)];
    rerr = [rerr,abs(ave(ind)+3)/3];
    plotyy(1:ind,ave,1:ind,rerr); drawnow;
    disp([ind,ave(ind)])
end
% hist(E)

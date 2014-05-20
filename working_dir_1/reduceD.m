function [mpsB] = reduceD(mpsA,mpoX,DB,factorization1,factorization2)
%
% factorization1 is a string specifying the factorization of the initial
% guess; it can be {'svd',,'qr','fastsvd','fastqr'} 
% factorization2 is a string specifying the factorization during left and
% right sweeps; it can be {'svd','qr'} 
%
% [mpsB]=reduceD(mpsA,mpoX,DB,factorization)
% [mpsB,tempo,tempo_loop]=reduceD(mpsA,mpoX,DB,factorization,tempo,tempo_loop)

N = length(mpsA); 
mpsB = fast_truncate(mpomps(mpoX,mpsA),'rl',DB,factorization1); 

% initialization of the storage 
Cstorage = initCstorage(mpsB,mpoX,mpsA,N);

% ****************** cycle 1: j -> j+1 (from 1 to N-1) ****************
for j=1:(N-1)
    % optimization
    Cleft=Cstorage{j};
    Cright=Cstorage{j+1};
    A=mpsA{j}; X=mpoX{j}; 
    [B,~]=reduceD2_onesite(A,X,Cleft,Cright); 
    [B,U]=prepare_onesite(B,'lr',factorization2);
    mpsB{j}=B; 

    % storage-update
    Cstorage{j+1}=updateCleft(Cleft,B,X,A); 
end

% ****************** cycle 2: j -> j-1 (from N to 2) ******************
for j=N:(-1):2
    % optimization
    Cleft=Cstorage{j};
    Cright=Cstorage{j+1};
    A=mpsA{j}; X=mpoX{j}; 
    [B,~]=reduceD2_onesite(A,X,Cleft,Cright); 
    [B,U]=prepare_onesite(B,'rl',factorization2);
    mpsB{j}=B; 
    % storage-update
    Cstorage{j}=updateCright(Cright,B,X,A); 
end

mpsB{1}=contracttensors(mpsB{1},3,2,U,2,1); 
mpsB{1}=permute(mpsB{1},[1,3,2]);


    
    
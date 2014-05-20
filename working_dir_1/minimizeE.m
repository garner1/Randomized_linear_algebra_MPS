function [E,mps]=minimizeE(hset,D,precision,mpsB,factorization)
%
% factorization: 'svd' or 'qr'
%
[M,N]=size(hset); 
d=size(hset{1,1},1); 
mps=createrandommps(N,D,d); 
mps=prepare(mps,factorization);

% storage-initialization
Hstorage=initHstorage(mps,hset,d);
if ~isempty(mpsB), Cstorage=initCstorage(mps,[],mpsB,N); end
P=[];

% optimization sweeps 
while 1
    Evalues=[];

    % ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
    for j=1:(N-1)
        % projector-calculation 
        if ~isempty(mpsB)
            B=mpsB{j};
            Cleft=Cstorage{j};
            Cright=Cstorage{j+1}; 
            P=calcprojector_onesite(B,Cleft,Cright);
        end

        % optimization
        Hleft=Hstorage(:,j);
        Hright=Hstorage(:,j+1);
        hsetj=hset(:,j);
        [A,E]=minimizeE_onesite(hsetj,Hleft,Hright,P); 
        [A,U]=prepare_onesite(A,'lr',factorization);
        mps{j}=A; 
        Evalues=[Evalues,E];

        % storage-update 
        for m=1:M
            h=reshape(hset{m,j},[1,1,d,d]);
            Hstorage{m,j+1}=updateCleft(Hleft{m},A,h,A); 
        end
        if ~isempty(mpsB) 
            Cstorage{j+1}=updateCleft(Cleft,A,[],B);
        end
    end
    
    % ****************** cycle 2: j -> j-1 (from N to 2) ****************** 
    for j=N:(-1):2
        % projector-calculation 
        if ~isempty(mpsB)
            B=mpsB{j};
            Cleft=Cstorage{j};
            Cright=Cstorage{j+1}; 
            P=calcprojector_onesite(B,Cleft,Cright);
        end
        
        % minimization
        Hleft=Hstorage(:,j);
        Hright=Hstorage(:,j+1);
        hsetj=hset(:,j); 
        [A,E]=minimizeE_onesite(hsetj,Hleft,Hright,P); 
        [A,U]=prepare_onesite(A,'rl',factorization);
        mps{j}=A; 
        Evalues=[Evalues,E];

        % storage-update 
        for m=1:M
            h=reshape(hset{m,j},[1,1,d,d]);
            Hstorage{m,j}=updateCright(Hright{m},A,h,A); 
        end
        if ~isempty(mpsB) 
            Cstorage{j}=updateCright(Cright,A,[],B);
        end
    end
    if (std(Evalues)/abs(mean(Evalues))<precision) 
        mps{1}=contracttensors(mps{1},3,2,U,2,1); 
        mps{1}=permute(mps{1},[1,3,2]);
        break;
    end
end






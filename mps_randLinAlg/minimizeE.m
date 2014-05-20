function [E,mps,Elist]=minimizeE(hset,D,precision,mpsB,its,l,prob)

[M,N]=size(hset); 
d=size(hset{1,1},1); 
mps=createrandommps(N,D,d); 
%
% bring the state into right canonical form
%
mps=prepare(mps);

% storage-initialization
Hstorage=initHstorage(mps,hset,d);
if ~isempty(mpsB), Cstorage=initCstorage(mps,[],mpsB,N); end
P=[];

% optimization sweeps 
Elist = [];
while 1
% for ind = 1:10
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
        [A,E]=minimizeE_onesite(hsetj,Hleft,Hright,P,its,l);
        if isempty(Elist),
% 
%             TO PREPARE THE TENSOR YOU CAN USE QR OR SFFT
% 
            [A,U] = prepare_onesite(A,'lr');
            mps{j} = A; 
            Evalues = [Evalues,E];
            Elist = [Elist,E];
        end
        Eold = Elist(end);
        if E < Eold,
            [A,U] = prepare_onesite(A,'lr');
            mps{j} = A; 
            Evalues = [Evalues,E];
            Elist = [Elist,E];
        end
        if E >= Eold, 
            if rand() < prob, 
                [A,U] = prepare_onesite(A,'lr');
                mps{j} = A; 
                Evalues = [Evalues,E];
                Elist = [Elist,E];
            else
                Evalues=[Evalues,Eold];
                Elist = [Elist,Eold];
                [A,U] = prepare_onesite(mps{j},'lr');
            end
        end
plot(Elist);
drawnow;
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
        [A,E]=minimizeE_onesite(hsetj,Hleft,Hright,P,its,l);
        Eold = Evalues(end);
        if E < Eold,
            [A,U] = prepare_onesite(A,'rl');
            mps{j} = A; 
            Evalues = [Evalues,E];
            Elist = [Elist,E];
        end
        if E >= Eold, 
            if rand() < prob, 
                [A,U] = prepare_onesite(A,'rl');
                mps{j} = A; 
                Evalues = [Evalues,E];
                Elist = [Elist,E];
            else
                Evalues = [Evalues,Eold];
                Elist = [Elist,Eold];
                [A,U] = prepare_onesite(mps{j},'rl');
            end
        end
plot(Elist);
drawnow;
        % storage-update 
        for m = 1:M
            h = reshape(hset{m,j},[1,1,d,d]);
            Hstorage{m,j} = updateCright(Hright{m},A,h,A); 
        end
        if ~isempty(mpsB) 
            Cstorage{j} = updateCright(Cright,A,[],B);
        end
    end
    if (std(Evalues)/abs(mean(Evalues))<precision) 
        mps{1} = contracttensors(mps{1},3,2,U,2,1); 
        mps{1} = permute(mps{1},[1,3,2]);
        E = Elist(end);
        break;
    end
end






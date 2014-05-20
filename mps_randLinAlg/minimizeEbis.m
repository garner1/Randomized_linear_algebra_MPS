function [E,mps]=minimizeEbis(hset,D,precision,mpsB)

[M,N]=size(hset); 
d=size(hset{1,1},1); 
mps=createrandommps(N,D,d); 
mps=prepare(mps);

% storage-initialization
% this is time consuming; find a way to speed it up
Hstorage=initHstorage(mps,hset,d);

% optimization sweeps 
while 1
    Evalues=[];
    % ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
    for j=1:(N-1)
        % optimization
        Hleft=Hstorage(:,j);
        Hright=Hstorage(:,j+1);
        hsetj=hset(:,j); 
        [A,E]=minimizeE_onesite(hsetj,Hleft,Hright,[]); 
        [A,U]=prepare_onesite(A,'lr');
        mps{j}=A; 
        Evalues=[Evalues,E];
        mps{j+1}=contracttensors(mps{j+1},3,1,U,2,2); 
        mps{j+1}=permute(mps{j+1},[3,1,2]);

        % storage-update 
        for m=1:M
            h=reshape(hset{m,j},[1,1,d,d]);
            Hstorage{m,j+1}=updateCleft(Hleft{m},A,h,A); 
        end
% disp(abs(real(E)))
%         if abs(real(E))<threshold, break;end
    end
%     if abs(real(E))<threshold,break;end
    
    % ****************** cycle 2: j -> j-1 (from N to 2) ****************** 
    for j=N:(-1):2
        % minimization
        Hleft=Hstorage(:,j);
        Hright=Hstorage(:,j+1);
        hsetj=hset(:,j); 
        [A,E]=minimizeE_onesite(hsetj,Hleft,Hright,[]); 
        [A,U]=prepare_onesite(A,'rl');
        mps{j}=A; 
        Evalues=[Evalues,E];
        mps{j-1}=contracttensors(mps{j-1},3,2,U,2,1); 
        mps{j-1}=permute(mps{j-1},[1,3,2]);
        % storage-update 
        for m=1:M
            h=reshape(hset{m,j},[1,1,d,d]);
            Hstorage{m,j}=updateCright(Hright{m},A,h,A); 
        end
% disp(abs(real(E)))
disp(std(Evalues)/abs(mean(Evalues)))
%         if abs(real(E))<threshold,break;end
    end
%     if abs(real(E))<threshold,break;end

end






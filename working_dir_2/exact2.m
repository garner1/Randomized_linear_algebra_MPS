function [p,pp,u,s,v] = exact2(mpox,mps,site)       
%
%this is used to check consistency with tensor network approach
%
% mpox=mpoid;
if site==1,
    N = size(mpox,2);
    [x,~,~,~] = size(mpox{1});
    [m,~,~] = size(mps{1});

    T = contracttensors(mpox{N},4,4,mps{N},3,3);
    T = contracttensors(T,5,[3],conj(mps{N}),3,[3]);
    T0 = permute(T,[5,1,3,6,2,4]);

    for ind = N-1:-1:2    
        T = contracttensors(mpox{ind},4,4,mps{ind},3,3);
        T = contracttensors(T,5,[3],conj(mps{ind}),3,[3]);
        T = permute(T,[5,1,3,6,2,4]);
        T0 = contracttensors(T,6,[4,5,6],T0,6,[1,2,3]);
    end
    p = reshape(T0,[m*x*m,m*x*m]);
    pp = reshape(permute(reshape(p,[m,m,m,m]),[1,3,2,4]),[m*m,m*m]);
    [u,s,v]=svd(p);
else
    N = size(mpox,2);
    [x,~,~,~] = size(mpox{1});
    [m,~,~] = size(mps{1});
    T = contracttensors(mpox{site-1},4,4,mps{site-1},3,3);
    T = contracttensors(T,5,[3],conj(mps{site-1}),3,[3]);
    T0 = permute(T,[5,1,3,6,2,4]);
    for ind = site-2:-1:1    
        T = contracttensors(mpox{ind},4,4,mps{ind},3,3);
        T = contracttensors(T,5,[3],conj(mps{ind}),3,[3]);
        T = permute(T,[5,1,3,6,2,4]);
        T0 = contracttensors(T,6,[4,5,6],T0,6,[1,2,3]);
    end
    for ind = N:-1:site+1    
        T = contracttensors(mpox{ind},4,4,mps{ind},3,3);
        T = contracttensors(T,5,[3],conj(mps{ind}),3,[3]);
        T = permute(T,[5,1,3,6,2,4]);
        T0 = contracttensors(T,6,[4,5,6],T0,6,[1,2,3]);
    end
    p = reshape(T0,[m*x*m,m*x*m]);
    pp = reshape(permute(reshape(p,[m,m,m,m]),[1,3,2,4]),[m*m,m*m]);
    [u,s,v]=svd(p);
end    


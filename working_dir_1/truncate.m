function [mps,U]=truncate(mps,direction,maxD,factorization) 

% this function uses svd to truncate the mps
% since the singular values are normalized it conserves the norm of the
% input state
N=length(mps);
switch direction    
    case 'lr'
        for i=1:N-1 
            [mps{i},U]=truncate_onesite(mps{i},'lr',maxD,factorization);
            mps{i+1}=contracttensors(U,2,2,mps{i+1},3,1);
        end
        [mps{N},U]=truncate_onesite(mps{N},'lr',maxD);
        
    case 'rl'
        for i=N:-1:2 
            [mps{i},U]=truncate_onesite(mps{i},'rl',maxD,factorization);
            mps{i-1}=contracttensors(mps{i-1},3,2,U,2,1);
            mps{i-1}=permute(mps{i-1},[1,3,2]); 
        end
        [mps{1},U]=truncate_onesite(mps{1},'rl',maxD,factorization);
end


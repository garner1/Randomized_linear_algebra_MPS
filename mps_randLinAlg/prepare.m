function [mps]=prepare(mps) 
%preprare a R-canonical state
N=length(mps);


for i=N:-1:2 
    [mps{i},U]=prepare_onesite(mps{i},'rl');
    mps{i-1}=contracttensors(mps{i-1},3,2,U,2,1);
    mps{i-1}=permute(mps{i-1},[1,3,2]); 
end


function [mps2] = mpomps(mpo,mps)

N=length(mps);
mps2=cell(1,N);

for indN=1:N
    [mpoD1,mpoD2,~,d]=size(mpo{indN});
    [mpsD1,mpsD2,~]=size(mps{indN});
    
    %[mpoD1,mpoD2,d,d][mpsD1,mpsD2,d]->[mpoD1,mpoD2,d,mpsD1,mpsD2]
    mps2{indN}=contracttensors(mpo{indN},4,4,mps{indN},3,3);

    %[mpoD1,mpoD2,d,mpsD1,mpsD2]->[mpoD1,mpsD1,mpoD2,mpsD2,d]
    mps2{indN}=permute(mps2{indN},[1 4 2 5 3]);
    mps2{indN}=reshape(mps2{indN},[mpoD1*mpsD1,mpoD2*mpsD2,d]);
   
end


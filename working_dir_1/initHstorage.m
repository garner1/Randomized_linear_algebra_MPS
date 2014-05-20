function [Hstorage]=initHstorage(mps,hset,d)

[M,N]=size(hset);
Hstorage=cell(M,N+1);
for m=1:M, Hstorage{m,1}=1; Hstorage{m,N+1}=1; end 
for j=N:-1:2
    for m=1:M
        h=reshape(hset{m,j},[1,1,d,d]); 
        Hstorage{m,j}=updateCright(Hstorage{m,j+1},mps{j},h,mps{j});
    end
end






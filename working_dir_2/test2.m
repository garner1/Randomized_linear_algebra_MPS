clear all;
tlist = [];
for n = 1000:1000:10000
    A = randn(n);
    tic
    B = reshape(A,[n^2,1]);
    tlist = [tlist,toc];
end

plot(tlist)
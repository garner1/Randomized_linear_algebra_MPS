N = 50;
D = 2^8;
d = 2;

mps00=createrandommps(N,D,d); 
its = [];
tic
mps0=prepare(mps00,its);
toc

its = 5;
tic
mps1=prepare(mps00,its);
toc

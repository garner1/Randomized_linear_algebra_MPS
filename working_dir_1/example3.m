N=20;D=5;d=2; 

mpsA=createrandommps(N,D,d);
mpsB=myrndmps(N,D,d,'right');

oset=cell(1,N);
id=eye(d); 
for j=1:N, oset{1,j}=id; end;

mz=expectationvalue(mpsB,oset);
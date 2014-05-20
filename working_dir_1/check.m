clear all;
N = 20;
d = 2;
D = 2^7;
mpo = mpo_hh(N,1,1,1,0,0);
mpsA = randmps_norm(N,D);
mpsA1 = reduceD(mpsA,mpo,D,'fastqr','qr');
mpsA2 = reduceD(mpsA,mpo,D,'svd','svd');

norm1 = overlap(mpsA1,mpsA1);
norm2 = overlap(mpsA2,mpsA2);

disp(real(overlap(mpsA1,mpsA2)/sqrt(norm1*norm2)))
% % 
% %check that the state is in rigth canonical form
% %
% ind = 2;
% mpsA1{ind}(:,:,1)*mpsA1{ind}(:,:,1)'+mpsA1{ind}(:,:,2)*mpsA1{ind}(:,:,2)'
% mpsA2{ind}(:,:,1)*mpsA2{ind}(:,:,1)'+mpsA2{ind}(:,:,2)*mpsA2{ind}(:,:,2)'

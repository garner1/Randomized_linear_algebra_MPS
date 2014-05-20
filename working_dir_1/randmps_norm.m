function [mps]=randmps_norm(N,D)
%PREPARE A LEFT NORMALIZED INHOMOGENEOUS RANDOM MPS STATE. The output is called mps. THE
%INPUT ARGUMENTS ARE THE NUMBER OF QUBITS AND THE MAX AUX DIMENSION.

%TO HAVE A NORMALIZED STATE ONE HAS TO USE DIFFERENT DIMENSIONS FOR THE 
% A MATRICES: THE FIRST LEFT HAS TO BE dxd, the second dxD1, then D1xD2,...
% D(N-1)xD(N)
% 1xd1x2 d1xd2x2 d2x1x2
% 2xd1 2d1xd2 2d2x1
% U(d1=2,d1) U(d2=4,d2) U(d2,1)x2
% for left normalization each U has to be an isometry such that U'*U=id, 
% where U can be obtained from a unitary operator truncating some of the
% columns

mps=cell(1,N); % each entry is tensor with 3 indices: first two are the indices of the A matrix (auxiliary space). The last is
%is the physical space

for ind=1:ceil(log2(D))-1
    U=randU(2^ind);
    U=reshape(U,[2,2^(ind-1),2^ind]);
    mps{ind}=permute(U,[2,3,1]);
end

U=randU(2^ceil(log2(D)));
U=U(:,1:D);
U=reshape(U,[2,2^(ceil(log2(D))-1),D]);
mps{ceil(log2(D))}=permute(U,[2,3,1]);

for ind=ceil(log2(D))+1:N-1
    U=randU(2*D);
    U=U(:,1:D);
    U=reshape(U,[2,D,D]);
    mps{ind}=permute(U,[2,3,1]);
end
U=randU(D)/sqrt(2); %the sqrt(2) is used because in the final contraction one has the sum of d=2 1's
mps{N}(:,1,1)=U(:,1);mps{N}(:,1,2)=U(:,2);

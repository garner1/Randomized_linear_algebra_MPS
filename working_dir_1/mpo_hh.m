function [mpo] = mpo_hh(N,jx,jy,jz,h,lambda)

% CREATE THE MPO OPERATOR FOR THE HEISENBERG HAMILTONIAN WITH N QUBITS AND
% GIVEN PARAMETERS;
% H=\SUM_sites JX*X*X+JY*Y*Y+JZ*Z*Z+H*SZ-LAMBDA*I

mpo=cell(1,N);
d=2;
sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2); 

%this part correspond to left side of Eq.183 in Schollwock review
aux=zeros(1,5,d,d);
for d_index2=1:d
    for d_index1=1:d
        aux(1,:,d_index1,d_index2)=[h*sz(d_index1,d_index2)-(lambda/N)*id(d_index1,d_index2), jx*sx(d_index1,d_index2),...
        jy*sy(d_index1,d_index2), jz*sz(d_index1,d_index2), id(d_index1,d_index2)];
    end
end
mpo{1}=aux;
aux2=aux;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this part correspond to right side of Eq.183 in Schollwock review
aux=zeros(5,1,d,d);
for d_index2=1:d
    for d_index1=1:d
        aux(:,1,d_index1,d_index2)=[id(d_index1,d_index2), sx(d_index1,d_index2), sy(d_index1,d_index2),...
        sz(d_index1,d_index2), h*sz(d_index1,d_index2)-(lambda/N)*id(d_index1,d_index2)];
    end
end
mpo{N}=aux;
aux3=aux;
%%%%%%%%%%%%%%%%%%%%%%%

%this part correspond to generic non-boundary term in Eq.182 of Schollwock
%review
for indN=2:N-1
    aux=zeros(5,5,d,d);
    aux(5,:,:,:)=aux2;
    aux(:,1,:,:)=aux3;
    mpo{indN}=aux;
end
    

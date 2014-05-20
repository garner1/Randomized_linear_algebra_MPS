function [mpo] = mpo_hh(N,jx,jy,jz,h,phi)
% [L,R,T,B]
% CREATE THE MPO OPERATOR FOR THE HEISENBERG HAMILTONIAN WITH N QUBITS AND
% GIVEN PARAMETERS;
% H=\SUM_sites JX*X*X+JY*Y*Y+JZ*Z*Z+H*SZ
%
% taken from "EFFICIENT MPS ALGORITHM FOR PERIODIC BOUNDARY CONDITIONS AND APPLICATIONS" 
%

mpo=cell(1,N);
d=2;
sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2); 

mpo{1} = zeros(5,5,d,d);
%
%for twisted boundary conditions
%
opex = cos(phi)*sx+1i*sin(phi)*sy;
opey = -1i*sin(phi)*sx+1i*cos(phi)*sy;
for d_index2 = 1:d
    for d_index1 = 1:d
        mpo{1}(1,:,d_index1,d_index2)=...
            [h*sz(d_index1,d_index2), jx*sx(d_index1,d_index2),...
        jy*sy(d_index1,d_index2), jz*sz(d_index1,d_index2), id(d_index1,d_index2)];
        mpo{1}(:,5,d_index1,d_index2)=...
            [id(d_index1,d_index2), opex(d_index1,d_index2), opey(d_index1,d_index2),...
        sz(d_index1,d_index2), 0];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mpo_aux = zeros(5,5,d,d);
for d_index2 = 1:d
    for d_index1 = 1:d
        mpo_aux(5,:,d_index1,d_index2)=...
            [h*sz(d_index1,d_index2), jx*sx(d_index1,d_index2),...
        jy*sy(d_index1,d_index2), jz*sz(d_index1,d_index2), id(d_index1,d_index2)];
        mpo_aux(:,1,d_index1,d_index2)=...
            [id(d_index1,d_index2), sx(d_index1,d_index2), sy(d_index1,d_index2),...
        sz(d_index1,d_index2), h*sz(d_index1,d_index2)];
    end
end

for ind=2:N, mpo{ind} = mpo_aux;end
